#include <errno.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <locale.h>
#include <signal.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <gsl/gsl_rng.h>

typedef struct{
    double	*x;
    double	*y;
    double	*z;
} vector_array;

const double	Lx = 100;	/* размер ячейки вдоль оси x (в единицах sigma) */
const double	Ly = 100;	/* размер ячейки вдоль оси y (в единицах sigma) */
const double	Lz = 100;	/* размер ячейки вдоль оси z (в единицах sigma) */
gsl_rng		*rng = NULL;	/* генератор случайных чисел из библиотеки gsl  */


/* функция возвращает случайное число в диапазоне 0..1 */
double get_random_value(void){
    return gsl_rng_uniform(rng);
}

/* потенциал Леннард-Джонса                 */
/* r2 -- квадрат расстояния между частицами */
double calc_potential(double r2){
    double r2i  = 1.0 / r2;
    double r6i = r2i * r2i * r2i;

    return 4.0 * r6i * (r6i - 1.0);
}

/* сила действующая на частицы взаимодействующие по потенциалу Леннард-Джонса */
/* здесь, для скорости, вычисляется неполная сила. Для получения полной силы  */
/* результат надо умножить на sqrt(r2)                                        */
/* r2 -- квадрат расстояния между частицами                                   */
double calc_force(double r2){
    double r2i  = 1.0 / r2;
    double r6i = r2i * r2i * r2i;

    return 24 * r2i * r6i * (2.0 * r6i - 1.0);
}

/* функция вычисляет кинетическую энергию системы */
double calc_kinetic_energy(vector_array *v, int count){
    double	energy;
    int		i;

    energy = 0;
    for(i = 0; i < count; i++){
	energy += v->x[i] * v->x[i] +
		  v->y[i] * v->y[i] +
		  v->z[i] * v->z[i];
    }
    return 0.5 * energy;
}

/* функция вычисляет ускорения всех частиц и возвращает */
/* потенциальную энергию системы                        */
double calc_accel_and_potential(vector_array *r, vector_array *a, int count){
    int		i, j;
    double	potential_energy;
    double	dx, dy, dz;
    double	f, r2;

    /* обнулить все ускорения */
    memset(a->x, 0, count * sizeof(double));
    memset(a->y, 0, count * sizeof(double));
    memset(a->z, 0, count * sizeof(double));

    /* посчитать ускорения и потенциал */
    potential_energy = 0;
    for(i = 0; i < count - 1; i++){
	for(j = i + 1; j < count; j++){
	    dx = r->x[i] - r->x[j];
	    dy = r->y[i] - r->y[j];
	    dz = r->z[i] - r->z[j];

	    r2 = dx * dx + dy * dy + dz * dz;
	    f = calc_force(r2);
	    potential_energy += calc_potential(r2);

	    a->x[i] += f * dx;
	    a->y[i] += f * dy;
	    a->z[i] += f * dz;

	    /* используем третий закон Ньютона */
	    a->x[j] -= f * dx;
	    a->y[j] -= f * dy;
	    a->z[j] -= f * dz;
	}
    }
    return potential_energy;
}

/* выполнить один шаг молекулярной динамики,                     */
/* найти новые значения координат, скоростей и ускорений частиц, */
/* а также кинетическую и потенциальную энергии                  */
void time_step(vector_array *r, vector_array *v, vector_array *a, int count,
		double Twall, int collision[], double dt,
		double *kinetic_energy, double *potential_energy){

    int		i;
    double	Vmax = sqrt(3.0 * Twall);

    /* очистить информацию о прошлых столкновениях */
    memset(collision, 0, count * sizeof(int));

    /* находим новые координаты */
    for(i = 0; i < count; i++){
	r->x[i] += (v->x[i] + 0.5 * a->x[i] * dt) * dt;
	r->y[i] += (v->y[i] + 0.5 * a->y[i] * dt) * dt;
	r->z[i] += (v->z[i] + 0.5 * a->z[i] * dt) * dt;

	/* частица налипает на стенку и отлетает с температурой стенки */
	if (r->x[i] <  0)  { r->x[i] = 0; collision[i] = 1; }
	if (r->y[i] <  0)  { r->y[i] = 0; collision[i] = 1; }
	if (r->z[i] <  0)  { r->z[i] = 0; collision[i] = 1; }

	/* частица налипает на стенку отлетает с температурой стенки */
	if (r->x[i] >= Lx) { r->x[i] = Lx; collision[i] = 1; }
	if (r->y[i] >= Ly) { r->y[i] = Ly; collision[i] = 1; }
	if (r->z[i] >= Lz) { r->z[i] = Lz; collision[i] = 1; }
    }

    /* находим новые скорости (часть 1) */
    for(i = 0; i < count; i++){
	if (collision[i] == 0){
	    v->x[i] += 0.5 * a->x[i] * dt;
	    v->y[i] += 0.5 * a->y[i] * dt;
	    v->z[i] += 0.5 * a->z[i] * dt;
	}
    }

    /* находим новые ускорения и вычисляем новый потенциал */
    *potential_energy = calc_accel_and_potential(r, a, count);

    /* находим новые скорости (часть 2) */
    for(i = 0; i < count; i++){
	if (collision[i] == 0){
	    v->x[i] += 0.5 * a->x[i] * dt;
	    v->y[i] += 0.5 * a->y[i] * dt;
	    v->z[i] += 0.5 * a->z[i] * dt;
	}else{
	    /* произошло столкновение, вычислим новую скорость частицы */
	    double	phi   = 2 * M_PI * get_random_value();
	    double	theta = M_PI * get_random_value();

	    v->x[i] = Vmax * cos(theta);
	    v->y[i] = Vmax * sin(theta) * cos(phi);
	    v->z[i] = Vmax * sin(theta) * sin(phi);

	    /* частица должна лететь во внутрь обьема */
	    if ((r->x[i] <= 0.0) && (v->x[i] <= 0.0)) v->x[i] = - v->x[i];
	    if ((r->y[i] <= 0.0) && (v->y[i] <= 0.0)) v->y[i] = - v->y[i];
	    if ((r->z[i] <= 0.0) && (v->z[i] <= 0.0)) v->z[i] = - v->z[i];

	    if ((r->x[i] >= Lx) && (v->x[i] >= 0.0)) v->x[i] = - v->x[i];
	    if ((r->y[i] >= Ly) && (v->y[i] >= 0.0)) v->y[i] = - v->y[i];
	    if ((r->z[i] >= Lz) && (v->z[i] >= 0.0)) v->z[i] = - v->z[i];

	}
    }

    /* вычисляем новое значение кинетической энергии */
    *kinetic_energy = calc_kinetic_energy(v, count);
}

/* Функция накапливает данные для построения гистограммы распределения    */
/* компоменты скорости частиц. Данные накапливаются в массиве data длиной */
/* bar_number. Значения выходящие за пределы v_min и v_max будут отнесены */
/* к крайним столбцам гистограммы. При первом обращении массив data       */
/* должен быть инициализирован нулями.                                    */
void collect_vi_data(double v[], int count, double v_min, double v_max,
			int data[], int bar_number){
    int		i, j;
    double	coeff;

    /* вычислить значения столбцов гистограммы */
    coeff = bar_number / (v_max - v_min);
    for(i = 0; i < count; i++){
	j = floor((v[i] - v_min) * coeff);

	if (j < 0 ) j = 0;
	if (j >= bar_number) j = bar_number - 1;
	data[j]++;
    }
}

/* Функция накапливает данные для построения гистограммы распределения */
/* модуля скорости частиц. Данные накапливаются в массиве data длиной  */
/* bar_number. Значения выходящие за пределы v_max будут отнесены      */
/* к крайним столбцам гистограммы. При первом обращении массив data    */
/* должен быть инициализирован нулями.                                 */
void collect_v_data(vector_array *v, int count, double v_max,
			int data[], int bar_number){
    int		i, j;
    double	v2, coeff;

    /* вычислить значения столбцов гистограммы */
    coeff = bar_number / v_max;
    for(i = 0; i < count; i++){
	v2 = sqrt( v->x[i] * v->x[i] +
		   v->y[i] * v->y[i] +
		   v->z[i] * v->z[i] );

	j = floor(v2 * coeff);
	if (j < 0 ) j = 0;
	if (j >= bar_number) j = bar_number - 1;

	data[j]++;
    }
}

/* Сохранить в файл гистограмму распределения скоростей частиц. Имя файла  */
/* хранится в name, данные в data[], bar_number задает количество столбцов */
/* в гистограмме, mode -- режим использованный для построения гистограммы  */
/* (может принимать значения "vx", "vy", "vz" или "v")                     */
int save_speed_distribution(char *name, int data[], int bar_number,
		char *mode, double v_min, double v_max,
		int count, double kinetic_energy, double t){

    int		i;
    double	coeff;
    FILE	*f;

    coeff = (v_max - v_min) / bar_number;

    /* сохранить гистограмму в файл с именем из name */
    if ((f = fopen(name, "w")) == NULL){
	printf("Ошибка : %s\n", strerror(errno));
	return 0;
    }

    fprintf(f, "# mode=%s\n", mode);
    fprintf(f, "# count=%d\n", count);
    fprintf(f, "# Ek=%lf\n", kinetic_energy);
    fprintf(f, "# t=%lf\n", t);
    fprintf(f, "# ------------------------------------------\n");
    fprintf(f, "# v\tnumber\n");
    fprintf(f, "# ------------------------------------------\n");
    for(i = 0; i < bar_number; i++){
	fprintf(f, "%lf\t%d\n", v_min + (i + 0.5) * coeff , data[i]);
    }

    fclose(f);
    printf("Сохранено\n");
    return 1;
}

/* сохранить в файл координаты и скорости всех частиц, а также */
/* другие данные которые могут понадобиться для дальнейшего    */
/* продолжения вычислений. Имя файла хранится в name           */
int save_data_file(char *name,
		vector_array *r, vector_array *v, int count,
		double kinetic_energy, double potential_energy,
		double Twall, double t){

    int		i;
    FILE	*f;

    if ((f = fopen(name, "w")) == NULL){
	printf("Ошибка : %s\n", strerror(errno));
	return 0;
    }

    fprintf(f, "# t=%lf\n", t);
    fprintf(f, "# Ek=%lf\n", kinetic_energy);
    fprintf(f, "# Ep=%lf\n", potential_energy);
    fprintf(f, "# E=%lf\n",  kinetic_energy + potential_energy);
    fprintf(f, "# count=%d\n", count);
    fprintf(f, "# Twall=%lf\n", Twall);
    fprintf(f, "# ------------------------------------------\n");
    fprintf(f, "# rx\try\trz\tvx\tvy\tvz\n");
    fprintf(f, "# ------------------------------------------\n");
    for(i = 0; i < count; i++){
	fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
	    r->x[i], r->y[i], r->z[i], v->x[i], v->y[i], v->z[i]);
    }
    fclose(f);
    printf("Сохранено\n");
    return 1;
}

/* загрузить из файла данные для дальнейшего продолжения вычислений. */
/* Имя файла передается в name                                       */
int load_data_file(const char *name,
		int *count, vector_array *r, vector_array *v,
		double *Twall, double *t){

    int		i;
    char	s[80];
    FILE	*file;

    /* открыть файл */
    if ((file = fopen(name, "r")) == NULL){
	printf("Ошибка чтения файла %s : %s\n", name, strerror(errno));
	return 0;
    }

    /* прочитать время */
    if (fscanf(file, "# t=%lf\n", t) != 1) goto error;
    if (fscanf(file, "# Ek=%[^\n]\n", s) != 1) goto error;
    if (fscanf(file, "# Ep=%[^\n]\n", s) != 1) goto error;
    if (fscanf(file, "# E=%[^\n]\n",  s) != 1) goto error;

    /* прочитать количество частиц */
    if (fscanf(file, "# count=%d\n", count) != 1) goto error;
    /* прочитать температуру стенки */
    if (fscanf(file, "# Twall=%lf\n", Twall) != 1) goto error;
    if (fscanf(file, "# --%[^\n]\n",  s) != 1) goto error;
    if (fscanf(file, "# rx%[^\n]\n",  s) != 1) goto error;
    if (fscanf(file, "# --%[^\n]\n",  s) != 1) goto error;

    /* проверить правильность прочитанных данных */
    if ((*t < 0) || (*count <= 0) || (*Twall <= 0)) goto error;

    /* выделить память под координаты молекул */
    r->x = malloc(*count * sizeof(double));
    r->y = malloc(*count * sizeof(double));
    r->z = malloc(*count * sizeof(double));

    /* выделить память под скорости молекул */
    v->x = malloc(*count * sizeof(double));
    v->y = malloc(*count * sizeof(double));
    v->z = malloc(*count * sizeof(double));

    if ((r->x == NULL) || (r->y == NULL) || (r->z == NULL) ||
	(v->x == NULL) || (v->y == NULL) || (v->z == NULL)){
	fprintf(stderr, "Не хватает памяти\n");
	fclose(file);
	return 0;
    }

    /* прочитать координаты и скорости из файла */
    for(i = 0; i < *count; i++){
	int result = fscanf(file, "%lf%lf%lf%lf%lf%lf\n",
	    &r->x[i], &r->y[i], &r->z[i], &v->x[i], &v->y[i], &v->z[i]);
	if (result != 6) goto error;
    }
    fclose(file);
    printf("Успешно прочитали файл %s,  count=%d\n", name, *count);
    return 1;

  error:
    fclose(file);
    printf("Ошибка чтения файла %s\n", name);
    return 0;
}

/* функция задает начальные координаты и скорости частиц */
void initialize(vector_array *r, vector_array *v, int count, double Vmax){
    int		i;
    double	px, py, pz;
    int		nx, ny, nz;

    /* раскидаем частицы равномерно по всему обьему, */
    /* но сперва подберем размер ячейки              */
    nx = ceil(Lx / pow(count, 1.0 / 3.0));
    ny = ceil(Ly / pow(count, 1.0 / 3.0));
    nz = ceil(Lz / pow(count, 1.0 / 3.0));
    while (nx * ny * nz < count){
	nx++; ny++; nz++;
    }

    r->x[0] = 0.5 * (Lx / nx);
    r->y[0] = 0.5 * (Ly / ny);
    r->z[0] = 0.5 * (Lz / nz);

    /* складываем кубики в ящик */
    for(i = 1; i < count; i++){
	r->x[i] = r->x[i - 1] + (Lx / nx);
	r->y[i] = r->y[i - 1];
	r->z[i] = r->z[i - 1];
	if (r->x[i] > Lx){
	    r->x[i] = 0.5 * (Lx / nx);
	    r->y[i] += (Ly / ny);
	    if (r->y[i] > Ly){
		r->y[i] = 0.5 * (Ly / ny);
		r->z[i] += (Lz / nz);
	    }
	}
    }

    /* скорости частиц пускай будут случайными */
    px = 0; py = 0; pz = 0;
    for(i = 0; i < count; i++){
	px += v->x[i] = Vmax * (2 * get_random_value() - 1.0);
	py += v->y[i] = Vmax * (2 * get_random_value() - 1.0);
	pz += v->z[i] = Vmax * (2 * get_random_value() - 1.0);
    }

    /* суммарный импульс должен равняться нулю */
    for(i = 0; i < count; i++){
	v->x[i] -= px / count;
	v->y[i] -= py / count;
	v->z[i] -= pz / count;
    }
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

/* функция спрашивает у пользователя какое распределение скоростей     */
/* он желает получить, спрашивает количество столбцов в гистограмме    */
/* и куда сохранить результат, а также (опционально) рисует полученную */
/* гистограмму                                                         */
void speed_distribution(vector_array *v, int count,
		double kinetic_energy, double t){

    int		*data, ch, bar_number;
    double	*v_proection, v_min, v_max;
    char	*mode, name[80];

    while(1){
	printf("Выберите вид распределения\n");
	printf("  0. Вернуться назад\n");
	printf("  1. Распределение скоростей по оси X\n");
	printf("  2. Распределение скоростей по оси Y\n");
	printf("  3. Распределение скоростей по оси Z\n");
	printf("  4. Распределение модуля скорости\n");
	printf("\n");
	printf("Ваш выбор : ");

	while((ch = getchar()) < ' ');
	switch(ch){
	    case '0':	/* Выход */
		return;

	    case '1':	/* Распределение скоростей по оси X */
		mode = "vx";
		v_proection = v->x;
		break;

	    case '2':	/* Распределение скоростей по оси Y */
		mode = "vy";
		v_proection = v->y;
		break;

	    case '3':	/* Распределение скоростей по оси Z */
		mode = "vz";
		v_proection = v->z;
		break;

	    case '4':	/* Распределение модуля скорости    */
		mode = "v";
		v_proection = NULL;
		break;

	    default:
		printf("Не понимаю, попробуем еще раз...\n");
		continue;
	}

	for(bar_number = 0; bar_number <= 0; ){
	    printf("Введите число столбцов гистограммы : ");
	    if (scanf("%d", &bar_number) == 1) break;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}

	memset(name, 0, sizeof(name));
	printf("Введите имя файла : ");
	scanf("%79s", name);

	/* выделить память для гистограммы и инициализировать ее нулями */
	if ((data = malloc(bar_number * sizeof(int))) == NULL){
	    printf("Ошибка : %s\n", strerror(errno));
	    return;
	}
	memset(data, 0, bar_number * sizeof(int));

	if (strcmp(mode, "v") != 0){
	    v_min = - 3.0 * sqrt(2.0 / 3.0 * kinetic_energy / count);
	    v_max = - v_min;
	    collect_vi_data(v_proection, count, v_min, v_max, data, bar_number);
	}else{
	    v_min = 0.0;
	    v_max = 3.0 * sqrt(2.0 * kinetic_energy / count);
	    collect_v_data(v, count, v_max, data, bar_number);
	}

	if (save_speed_distribution(name, data, bar_number,
		mode, v_min, v_max,
		count, kinetic_energy, t) == 1){

	    printf("Нарисовать график (y/n) ? ");
	    while((ch = getchar()) < ' ');
	    if ((ch == 'y') || (ch == 'Y')){
		char cmd[256];
		snprintf(cmd, sizeof(cmd), "plot \"%s\" using 1:2 with impulses lw 3; pause mouse; exit;\n", name);
		run_gnuplot(cmd);
	    }
	}
	free(data);
    }
}

/* функция спрашивает у пользователя какое распределение частиц */
/* в пространстве он желает получить, спрашивает куда сохранить */
/* результат и рисует полученное распределение                  */
void space_distribution(vector_array *r, vector_array *v, int count,
		double kinetic_energy, double potential_energy,
		double Twall, double t){

    int		ch, col1, col2;
    char	name[80], cmd[256];

    memset(name, 0, sizeof(name));
    while(1){
	printf("Выберите вид распределения\n");
	printf("  0. Вернуться назад\n");
	printf("  1. Проекция распределения частиц на XY\n");
	printf("  2. Проекция распределения частиц на YZ\n");
	printf("  3. Проекция распределения частиц на XZ\n");
	printf("\n");
	printf("Ваш выбор : ");

	while((ch = getchar()) < ' ');
	switch(ch){
	    case '0':	/* Выход */
		return;

	    case '1':	/* Проекция распределения частиц на XY */
		col1 = 1;
		col2 = 2;
		break;

	    case '2':	/* Проекция распределения частиц на YZ */
		col1 = 2;
		col2 = 3;
		break;

	    case '3':	/* Проекция распределения частиц на XZ */
		col1 = 1;
		col2 = 3;
		break;

	    default:
		printf("Не понимаю, попробуем еще раз...\n");
		continue;
	}

	if (*name == '\0'){
	    printf("Введите имя файла : ");
	    scanf("%79s", name);
	    if (save_data_file(name, r, v, count,
				kinetic_energy,
				potential_energy,
				Twall, t) != 1){
		memset(name, 0, sizeof(name));
		continue;
	    }
	}

	snprintf(cmd, sizeof(cmd),
	    "plot \"%s\" using %d:%d with points; pause mouse; exit;\n",
	    name, col1, col2);
	run_gnuplot(cmd);
    }
}

/* спросить имя файла и сохранить в файл координаты и скорости   */
/* всех частиц, а также другие данные которые могут понадобиться */
/* для дальнейшего продолжения вычислений.                       */
void save_snapshot(vector_array *r, vector_array *v, int count,
		double kinetic_energy, double potential_energy,
		double Twall, double t){

    char	name[80];

    memset(name, 0, sizeof(name));
    printf("Введите имя файла : ");
    scanf("%79s", name);

    save_data_file(name, r, v, count, kinetic_energy, potential_energy, Twall, t);
}

/* функция ввода нового значения шага по времени */
double input_new_dt(void){
    double	dt;

    for(dt = 0; dt <= 0; ){
	printf("Введите шаг по времени : ");
	if (scanf("%lf", &dt) == 1) break;
	while(getchar() != '\n');
	printf("Попробуйте еще раз\n");
    }
    return dt;
}

/* Для того чтобы попасть в эту фунцию надо нажать Ctrl-C в то время,     */
/* когда программа занята вычислениями. Функция спрашивает у пользователя */
/* чего он собственно хочет и выполняет его пожелания                     */
int action_request(vector_array *r, vector_array *v, int count,
		double kinetic_energy, double potential_energy,
		double Twall, double t, double *dt){
    int		ch;

    while(1){
	printf("\n");
	printf("Что будем делать, хозяин?\n");
	printf("  0. Продолжим вычисления\n");
	printf("  1. Посмотрим на распределение скоростей\n");
	printf("  2. Посмотрим на распределение частиц в пространстве\n");
	printf("  3. Сохраним данные в файл\n");
	printf("  4. Введем новый шаг по времени\n");
	printf("  5. Выход\n");
	printf("\n");
	printf("Ваш выбор : ");

	while((ch = getchar()) < ' ');
	switch(ch){
	    case '0':	/* Продолжим вычисления */
		return 0;

	    case '1':	/* Посмотрим на распределение скоростей */
		speed_distribution(v, count, kinetic_energy, t);
		continue;

	    case '2':	/* Посмотрим на распределение частиц в пространстве */
		space_distribution(r, v, count, kinetic_energy, potential_energy, Twall, t);
		continue;

	    case '3':	/* Сохраним данные в файл */
		save_snapshot(r, v, count, kinetic_energy, potential_energy, Twall, t);
		continue;

	    case '4':	/* Введем новый шаг по времени */
		printf("Текущее значение шага по времени dt=%lf\n", *dt);
		*dt = input_new_dt();
		continue;

	    case '5':	/* Выход */
		return 1;

	    default:
		printf("Не понимаю, попробуем еще раз...\n");
		continue;
	}
    }
    return 0;
}

/* функция запрещает/восстанавливает обычное поведение программы */
/* при нажатии на Ctrl-C                                         */
void set_ctrl_c_reaction(int enable){
    sigset_t		signal_set;

    sigemptyset(&signal_set);
    sigaddset(&signal_set, SIGINT);
    sigprocmask(enable ? SIG_UNBLOCK : SIG_BLOCK, &signal_set, NULL);
}

/* функция выводит подсказку о том как пользоваться этой программой */
void help(char *name){
    if (strrchr(name, '/') != NULL) name = strrchr(name, '/') + 1;
    printf("Use: %s [--load file]\n", name);
}

/* с этой функции начинается выполнение программы */
int main(int argc, char *argv[]){
    time_t		last_time;
    int			load_flag = 0;
    int			count;			/* количество частиц            */
    double		Twall;			/* температура стенки           */
    double		Vmax;			/* максимальная скорость частиц */
    double		t, dt;			/* время и шаг по времени       */
    vector_array	r;			/* координаты частиц            */
    vector_array	v;			/* скорости частиц              */
    vector_array	a;			/* ускорение частиц             */
    int			*collision;
    double		kinetic_energy;
    double		potential_energy;

    /* инициализируем генератор случайных чисел */
    rng = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(rng, time(NULL));

    /* напечатаем приветствие */
    printf("Моделирование канонического ансамбля\n\n");

    /* проанализируем аргументы командной строки */
    while(1){
	int			c;
	int			option_index	= 0;
	static struct option	long_options[]	= { { "help", 0, 0, 'h' },
						    { "load", 1, 0, 'l' },
						    { NULL,   0, 0, 0   } };

	c = getopt_long(argc, argv, "hl:", long_options, &option_index);
	if (c == -1) break;

	switch(c){
	    case 'h':	/* get help */
		help(argv[0]);
		exit(EXIT_SUCCESS);

	    case 'l':	/* load file */
		if ((optarg == NULL) || (*optarg == '\0')){
		    printf("Не хватает имени файла\n");
		    help(argv[0]);
		    exit(EXIT_FAILURE);
		}
		if (load_data_file(optarg, &count, &r, &v, &Twall, &t) != 1){
		    exit(EXIT_FAILURE);
		}
		load_flag = 1;
		break;

	    default:
		printf("Неизвестная опция : %c\n", c);
		help(argv[0]);
		exit(EXIT_FAILURE);
	};
    }

    if ( ! load_flag){
	for(count = 0; count <= 0; ){
	    printf("Введите число молекул : ");
	    if (scanf("%d", &count) == 1) break;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}

	for(Twall = 0; Twall <= 0; ){
	    printf("Введите температуру стенки : ");
	    if (scanf("%lf", &Twall) == 1) break;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}

	for(Vmax = 0; Vmax <= 0; ){
	    printf("Введите максимальную скорость частиц : ");
	    if (scanf("%lf", &Vmax) == 1) break;
	    while(getchar() != '\n');
	    printf("Попробуйте еще раз\n");
	}

	/* выделить память под координаты молекул */
	r.x = malloc(count * sizeof(double));
	r.y = malloc(count * sizeof(double));
	r.z = malloc(count * sizeof(double));

	/* выделить память под скорости молекул */
	v.x = malloc(count * sizeof(double));
	v.y = malloc(count * sizeof(double));
	v.z = malloc(count * sizeof(double));

	if ((r.x == NULL) || (r.y == NULL) || (r.z == NULL) ||
	    (v.x == NULL) || (v.y == NULL) || (v.z == NULL)){
	    fprintf(stderr, "Не хватает памяти, попробуйте уменьшить количество молекул\n");
	    exit(EXIT_FAILURE);
	}

	/* задать время, начальные координаты и скорости */
	t = 0;
	initialize(&r, &v, count, Vmax);
    }

    /* выделить память под ускорение молекул */
    a.x = malloc(count * sizeof(double));
    a.y = malloc(count * sizeof(double));
    a.z = malloc(count * sizeof(double));

    /* выделить память для запоминания факта столкновения со стенкой */
    collision = malloc(count * sizeof(int));


    if ((collision == NULL) ||
	(a.x == NULL) || (a.y == NULL) || (a.z == NULL)){

	fprintf(stderr, "Не хватает памяти, попробуйте уменьшить количество молекул\n");
	exit(EXIT_FAILURE);
    }

    /* посчитать ускорения и энергию системы */
    kinetic_energy   = calc_kinetic_energy(&v, count);
    potential_energy = calc_accel_and_potential(&r, &a, count);

    printf("-------------------------------------\n");
    printf("Среднее расстояние между частицами : %lf\n", pow(Lx * Ly * Lz / count, 1.0/3.0));
    printf("Начальная кинетическая энергия     : %lf\n", kinetic_energy);
    printf("Начальная потенциальная энергия    : %lf\n", potential_energy);
    printf("Время начала отсчёта               : %lf\n", t);
    printf("Температура стенки                 : %lf\n", Twall);
    printf("Полная энергия системы             : %lf\n", kinetic_energy + potential_energy);
    printf("-------------------------------------\n");

    /* вводим шаг по времени */
    dt = input_new_dt();

    last_time = time(NULL);
    set_ctrl_c_reaction(0);
    while(1){
	time_t		current_time = time(NULL);
	siginfo_t	siginfo;
	sigset_t	signal_set;
	struct timespec	timeout;

	t += dt;
	/* делаем один шаг метода молекулярной динамики */
	time_step(&r, &v, &a, count,
	    Twall, collision, dt, &kinetic_energy, &potential_energy);
	if (current_time != last_time){
	    last_time = current_time;
	    printf("\rt=%lf,  Ek=%lf,  Ep=%lf,  E=%lf,  T=%lf       ",
		t, kinetic_energy, potential_energy,
		kinetic_energy + potential_energy,
		(2.0 / 3.0) * kinetic_energy / count);
	    fflush(stdout);
	}

	/* ловим нажатие Ctrl-C */
	sigemptyset(&signal_set);
	sigaddset(&signal_set, SIGINT);
	timeout.tv_sec  = 0;
	timeout.tv_nsec = 0;
	if (sigtimedwait(&signal_set, &siginfo, &timeout) > 0){
	    /* нажали Ctrl-C */

	    printf("\rt=%lf,  Ek=%lf,  Ep=%lf,  E=%lf,  T=%lf       ",
		t, kinetic_energy, potential_energy,
		kinetic_energy + potential_energy,
		(2.0 / 3.0) * kinetic_energy / count);
	    fflush(stdout);

	    set_ctrl_c_reaction(1);
	    printf("\n");
	    /* выводим меню и спрашиваем чего хочет пользователь */
	    if (action_request(&r, &v, count,
		kinetic_energy, potential_energy,
		Twall, t, &dt) == 1){
		/* завершаем работу программы */
		break;
	    }
	    set_ctrl_c_reaction(0);
	}
    }
    printf("\n");

    return 0;
}