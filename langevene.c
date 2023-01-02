#include <errno.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <gsl/gsl_rng.h>

struct P{
    double	x;
    double	y;
    double	vx;
    double	vy;
};

gsl_rng		*rng = NULL;	/* генератор случайных чисел из библиотеки gsl  */


/* функция возвращает случайное число в диапазоне 0..1 */
double get_random_value(void){
    return gsl_rng_uniform(rng);
}

/* вычисляем траекторию броуновской частицы */
void langevene_solution(double F, int N, double dt, struct P *p0, struct P p[]){
    int		i;

    memcpy(&p[0], p0, sizeof(struct P));
    for(i = 1; i <= N; i++){
	p[i].x  = p[i - 1].x + p[i - 1].vx * dt;
	p[i].y  = p[i - 1].y + p[i - 1].vy * dt;
	p[i].vx = (1.0 - dt) * p[i - 1].vx + F * dt * (1.0 - 2.0 * get_random_value());
	p[i].vy = (1.0 - dt) * p[i - 1].vy + F * dt * (1.0 - 2.0 * get_random_value());
    }
}

/* сохранить в файл данные статистики */
void save_stat_data(FILE *f, int N, double dt, double F,
		int count, struct P pa[], double r2[], double v2[]){

    int		i;

    fprintf(f, "# N=%d\n", N);
    fprintf(f, "# dt=%lf\n", dt);
    fprintf(f, "# F=%lf\n", F);
    fprintf(f, "# count=%d\n", count);
    fprintf(f, "# ---------------------------------------------------\n");
    fprintf(f, "# t\t<x>\t<y>\t<vx>\t<vy>\t<x^2>\t<v^2>\n");
    fprintf(f, "# ---------------------------------------------------\n");
    for(i = 0; i <= N; i++)
	fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
	    i * dt, pa[i].x / count,  pa[i].y / count,
	            pa[i].vx / count, pa[i].vy / count,
	            r2[i] / count, v2[i] / count);
}

/* Запустить программу gnuplot для построения графиков.              */
/* Команда для программы gnuplot передается в cmd.                   */
/* Не забудьте про паузу в cmd, если хотите увидеть график на экране */
int run_gnuplot(const char *cmd){
    int		status;
    int		filedes[2];
    pid_t	pid;

    /* создаем канал передачи данных */
    if (pipe(filedes) < 0){
	printf("Ошибка : %s\n", strerror(errno));
	return 0;
    }

    /* делимся */
    if ((pid = fork()) == (pid_t) (-1)){
	printf("Ошибка : %s\n", strerror(errno));
	return 0;
    }

    /* деление прошло успешно */
    if (pid == 0){
	/* этот код будет выполняться ТОЛЬКО в дочернем процессе */

	/* перенаправляем стандатный ввод дочернего процесса на */
	/* вывод канала передачи данных                         */
	dup2(filedes[0], 0);

	/* закрываем ненужные более ресурсы */
	close(filedes[0]);
	close(filedes[1]);

	/* вызываем gnuplot */
	execlp("gnuplot", "gnuplot", NULL);

	/* если gnuplot запустился, мы до этого места уже не дойдем */
	/* ну а если мы сюда попали -- значит нужно умереть         */
	printf("Ошибка : %s\n", strerror(errno));
	exit(EXIT_FAILURE);
    }

    /* закрываем вывод канала передачи данных (нужен только в потомке) */
    close(filedes[0]);
    /* записываем команду в канал передачи данных, команда попадет на */
    /* стандартный вход программы gnuplot и будет исполнена им        */
    write(filedes[1], cmd, strlen(cmd));
    /* полностью закрываем канал передачи данных */
    close(filedes[1]);

    /* ждем окончания работы дочернего процесса */
    wait(&status);
    return status;
}

/* спросить имя файла и сохранить туда данные статистики */
void save_stat(int N, double dt, double F,
		int count, struct P pa[], double r2[], double v2[]){

    FILE	*f;
    char	name[80];

    memset(name, 0, sizeof(name));
    printf("Введите имя файла : ");
    scanf("%79s", name);

    if ((f = fopen(name, "w")) == NULL){
	printf("Не могу сохранить статистику в файл %s\n", name);
	return;
    }
    save_stat_data(f, N, dt, F, count, pa, r2, v2);
    fclose(f);
    printf("Сохранено\n");
}

/* Для того чтобы попасть в эту фунцию надо нажать Ctrl-C в то время,     */
/* когда программа занята вычислениями. Функция спрашивает у пользователя */
/* чего он собственно хочет и выполняет его пожелания                     */
int action_request(int N, double dt, double F, struct P p[],
		int count, struct P pa[], double r2[], double v2[]){

    char	stat_template[] = "/tmp/langevene-XXXXXX";
    char	pos_template[] = "/tmp/langevene-pos-XXXXXX";
    char	cmd[256];
    int		fd;
    FILE	*f;
    int		i, ch;

    /* инициализируем stat_template для создания временного файла */
    memset(stat_template + sizeof(stat_template) - 7, 'X', 6);
    stat_template[sizeof(stat_template) - 1] = '\0';

    /* создаем временный файл и пишем туда данные статистики */
    if ((fd = mkstemp(stat_template)) == -1){
	printf("Не могу создать временный файл для сохранения статистики\n");
	exit(EXIT_FAILURE);
    }
    if ((f = fdopen(fd, "w")) == NULL){
	close(fd);
	unlink(stat_template);
	printf("Не могу сохранить статистику во временный файл\n");
	exit(EXIT_FAILURE);
    }
    save_stat_data(f, N, dt, F, count, pa, r2, v2);
    fclose(f);

    /* инициализируем pos_template для создания временного файла */
    memset(pos_template + sizeof(pos_template) - 7, 'X', 6);
    pos_template[sizeof(pos_template) - 1] = '\0';

    /* создаем временный файл и пишем туда последнюю вычисленную траекторию */
    if ((fd = mkstemp(pos_template)) == -1){
	printf("Не могу создать временный файл для сохранения траектории\n");
	unlink(stat_template);
	exit(EXIT_FAILURE);
    }
    if ((f = fdopen(fd, "w")) == NULL){
	close(fd);
	unlink(pos_template);
	unlink(stat_template);
	printf("Не могу сохранить траекторию во временный файл\n");
	exit(EXIT_FAILURE);
    }
    for(i = 0; i <= N; i++) fprintf(f, "%lg\t%lg\t%lg\n", i * dt, p[i].x, p[i].y);
    fclose(f);

    while(1){
	printf("\n");
	printf("Что будем делать, хозяин?\n");
	printf("  1. Посмотрим на средний квадрат скорости частицы\n");
	printf("  2. Посмотрим на средний квадрат координаты частицы\n");
	printf("  3. Посмотрим на среднюю скорость частицы\n");
	printf("  4. Посмотрим на среднюю траекторию частицы\n");
	printf("  5. Посмотрим последнюю вычисленную траекторию частицы\n");
	printf("  6. Сохраним данные в файл\n");
	printf("  7. Выход\n");
	printf("\n");
	printf("Ваш выбор : ");

	while((ch = getchar()) < ' ');
	switch(ch){
	    case '1':	/* Посмотрим на средний квадрат скорости частицы */
		snprintf(cmd, sizeof(cmd), "plot '%s' using 1:7 with lines notitle; pause mouse; exit;\n", stat_template);
		run_gnuplot(cmd);
		continue;

	    case '2':	/* Посмотрим на средний квадрат координаты частицы */
		snprintf(cmd, sizeof(cmd), "plot '%s' using 1:6 with lines notitle; pause mouse; exit;\n", stat_template);
		run_gnuplot(cmd);
		continue;

	    case '3':	/* Посмотрим на среднюю скорость частицы */
		snprintf(cmd, sizeof(cmd), "plot '%s' using 4:5 with lines notitle; pause mouse; exit;\n", stat_template);
		run_gnuplot(cmd);
		continue;

	    case '4':	/* Посмотрим на среднюю траекторию частицы */
		snprintf(cmd, sizeof(cmd), "plot '%s' using 2:3 with lines notitle; pause mouse; exit;\n", stat_template);
		run_gnuplot(cmd);
		continue;

	    case '5':	/* Посмотрим последнюю вычисленную траекторию частицы */
		snprintf(cmd, sizeof(cmd), "plot '%s' using 2:3 with lines notitle; pause mouse; exit;\n", pos_template);
		run_gnuplot(cmd);
		continue;

	    case '6':	/* Сохраним данные в файл */
		save_stat(N, dt, F, count, pa, r2, v2);
		continue;

	    case '7':	/* Выход */
		unlink(stat_template);
		unlink(pos_template);
		return 1;

	    default:
		printf("Не понимаю, попробуем еще раз...\n");
		continue;
	}
    }
    return 0;
}

void print_status(int count, double epsilon){
    printf("\rcount=%d,  epsilon=%lf       ", count, epsilon);
    fflush(stdout);
}

/* функция выводит подсказку о том как пользоваться этой программой */
void help(char *name){
    if (strrchr(name, '/') != NULL) name = strrchr(name, '/') + 1;
    printf("Use: %s [--N count] [--dt dt] [--F max_force] [--eps epsilon]\n", name);
}

/* с этой функции начинается выполнение программы */
int main(int argc, char *argv[]){
    int			count;	/* количество итераций                  */
    int			N;	/* количество шагов по времени          */
    double		dt;	/* шаг по времени                       */
    double		F;	/* константа при случайной силе         */
    double		eps;	/* относительная точность               */
    struct P		p0;	/* начальные координаты и скорости      */
    struct P		*p;	/* массив координат и скоростей         */
    struct P		*pa;	/* массив средних координат и скоростей */
    double		*v2;	/* массив среднего квадрата скорости    */
    double		*r2;	/* массив среднего квадрата координаты  */

    /* инициализируем генератор случайных чисел */
    rng = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(rng, time(NULL));

    /* напечатаем приветствие */
    printf("Уравнение Ланжевена\n\n");

    N = 0;
    dt = 0;
    F = -1;
    eps = 0.05;

    /* проанализируем аргументы командной строки */
    while(1){
	int			c;
	int			option_index	= 0;
	static struct option	long_options[]	= { { "help", 0, 0, 'h' },
						    { "N",    1, 0, 'N' },
						    { "dt",   1, 0, 't' },
						    { "F",    1, 0, 'F' },
						    { "eps",  1, 0, 'e' },
						    { NULL,   0, 0, 0   } };

	c = getopt_long(argc, argv, "h", long_options, &option_index);
	if (c == -1) break;

	switch(c){
	    case 'h':	/* get help */
		help(argv[0]);
		exit(EXIT_SUCCESS);

	    case 'N':	/* количество шагов по времени */
		if ((optarg == NULL) || (*optarg == '\0')){
		    printf("Не хватает количества шагов по времени\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (sscanf(optarg, "%d", &N) != 1){
		    printf("Не могу прочитать количество шагов по времени\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (N <= 0){
		    printf("Количество шагов по времени слишком мало\n");
		    exit(EXIT_FAILURE);
		}
		break;

	    case 't':	/* шаг по времени */
		if ((optarg == NULL) || (*optarg == '\0')){
		    printf("Не хватает шага по времени\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (sscanf(optarg, "%lf", &dt) != 1){
		    printf("Не могу прочитать шаг по времени\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if ((dt <= 0) || (dt >= 1)){
		    printf("Шаг по времени должен быть меньше 1.0\n");
		    exit(EXIT_FAILURE);
		}
		break;

	    case 'F':	/* константа при случайной силе */
		if ((optarg == NULL) || (*optarg == '\0')){
		    printf("Не хватает константы при случайной силе\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (sscanf(optarg, "%lf", &F) != 1){
		    printf("Не могу прочитать константу при случайной силе\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (F < 0){
		    printf("Константа при случайной силе должна быть неотрицательна\n");
		    exit(EXIT_FAILURE);
		}
		break;

	    case 'e':	/* относительная точность */
		if ((optarg == NULL) || (*optarg == '\0')){
		    printf("Не хватает относительной точности\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (sscanf(optarg, "%lf", &eps) != 1){
		    printf("Не могу прочитать относительную точность\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if ((eps <= 0) || (eps > 0.5)){
		    printf("Относительная точность должна быть в интервале (0..0.5]\n");
		    exit(EXIT_FAILURE);
		}
		break;

	    default:
		printf("Неизвестная опция : %c\n", c);
		help(argv[0]);
		exit(EXIT_FAILURE);
	};
    }

    while(N <= 0){
	printf("Введите количество шагов по времени : ");
	if (scanf("%d", &N) != 1){
	    N = 0;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}
    }

    while((dt <= 0) || (dt >= 1.0)){
	printf("Введите шаг по времени : ");
	if (scanf("%lf", &dt) != 1){
	    dt = 0;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}
    }

    while(F < 0){
	printf("Введите константу при случайной силе : ");
	if (scanf("%lf", &F) != 1){
	    F = -1;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}
    }

    /* выделяем память */
    p  = calloc(N + 1, sizeof(struct P));
    pa = calloc(N + 1, sizeof(struct P));
    v2 = calloc(N + 1, sizeof(double));
    r2 = calloc(N + 1, sizeof(double));

    if ((p == NULL) || (pa == NULL) || (v2 == NULL) || (r2 == NULL)){
	printf("Не хватает памяти, попробуйте уменьшить количество шагов по времени\n");
	exit(EXIT_FAILURE);
    }

    /* зададим начальные координаты и скорости */
    p0.x = p0.y = p0.vx = p0.vy = 0.0;

    printf("-------------------------------------\n");
    printf("Количество шагов по времени        : %d\n", N);
    printf("Шаг по времени                     : %lf\n", dt);
    printf("Константа при случайной силе       : %lf\n", F);
    printf("-------------------------------------\n");

    count = 0;
    while(1){
	int		i;
	double		epsilon;

	/* вычисляем траекторию броуновской частицы */
	langevene_solution(F, N, dt, &p0, p);

	/* собираем статистику и вычисляем погрешность */
	epsilon = -1;
	for(i = 0; i <= N; i++){
	    double	rr2, vv2;
	
	    pa[i].x  += p[i].x;
	    pa[i].y  += p[i].y;
	    pa[i].vx += p[i].vx;
	    pa[i].vy += p[i].vy;
	    v2[i] += vv2 = p[i].vx * p[i].vx + p[i].vy * p[i].vy;
	    r2[i] += rr2 = p[i].x  * p[i].x  + p[i].y  * p[i].y;
	    count++;

	    if ((count > 1) && (v2[i] > 0) && (r2[i] > 0)){
		double ev = fabs((v2[i] - count * vv2) / ((count - 1) * v2[i]));
		double er = fabs((r2[i] - count * rr2) / ((count - 1) * r2[i]));
		if (ev > epsilon) epsilon = ev;
		if (er > epsilon) epsilon = er;
	    }
	}

	/* достигли требуемой точности? */
	if ((epsilon > 0) && (epsilon < eps)){
	    print_status(count, epsilon);
	    break;
	}

	if (count % 100 == 0) print_status(count, epsilon);
    }

    printf("\n");
    action_request(N, dt, F, p, count, pa, r2, v2);
    return 0;
}
