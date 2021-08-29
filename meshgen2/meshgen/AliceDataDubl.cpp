// В файле AliceDataDubl.cpp содержатся данные и функции
// уже реализованные в различных модулях программы AliceFlow.
// Поскольку сеточный генератор это самостооятельный проект то здесь необходимо продублировать данные и методы.
// Вбудущем при интеграцции сеточного генератора в программу AliceFlow данный файл включать ненадо.



#define Real double
#define Logical bool

// из файла AliceData.h
// begin

#define TEMP 0 // температура
#define VX 1 // горизонтальная скорость
#define VY 2 // вертикальная скорость
#define VAR_COUNT 3 // количество функций под которые требуется выделить память.

typedef struct TMYTASK_DATA {
	int maxnod; // максимальный номер узла (размерность массива)
    int maxelm; // максимально допустимое число элементов (элементов меньше чем узлов)
    // по умолчанию элементы треугольники. Тип элемента определяется при считывании файла с сеткой.
    int nve; // число узловых переменных элемента (программа работает только для квадратных).

    int **nvtx; // список узлов для каждого элемента
    // динамическое выделение памяти не используется
    Real *x, *y; // узловые координаты
    //int ***nvtxMG; // список узлов для каждого элемента для нескольких уровней вложенности.

    Logical **constr; // истина для фиксированных потенциалов (лог. переменная)
    Real **potent; // массив узловых потенциалов
    //Real *rthdsd; // правая часть системы уравнений

    // свойства материалов
    Real *rho, *cp, *lam;

	//Real **aelm; // для одного элемента A.
    //Real *rloc; // локальная правая часть

    int nodes; // число узлов в задаче
    int nelmts; // число элементов в модели

} MYTASK_DATA;

// end;

// файл myexporttecplot.cpp
// начало.

// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// универсально подходит и для треугольных и для квадратных ячеек.
void exporttecplotxy360(int nve, int nelmts, int nodes, int** nvtx, Real* x, Real* y, Real* potent)
{
	FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "tigel_interior_tec.dat"/*"temp_xy.PLT"*/, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"temperature\"\n");

		// запись имён переменных
		//fprintf(fp, "VARIABLES = x, y, \"temperature\" , Vx, Vy, Mag \n");
		if (btriangle) {
			fprintf(fp, "VARIABLES = x, y \n"); // без температуры (зачем лишние данные это же сеточный генератор).
		}
		else {
		   fprintf(fp, "VARIABLES = x, y, \"temperature\" \n");
		}

		// запись информации о зонах
		if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nodes, nelmts);
        if (nve==4) {
			if (btriangle) {
				// каждый четырёхугольник это два треугольника, поэтому 
				// треугольных элементов в два раза больше.
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nodes, 2*nelmts);
			}
			else {
				// истинные четырёхугольные элементы.
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", nodes, nelmts);
			}
		}

		int i=0; // счётчики 
		int j=0; // цикла for

		// запись x
	    for (i=0; i<nodes; i++) {	
				fprintf(fp, "%e ", x[i]);
				if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// запись y
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", y[i]);
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		if (!btriangle) {
		    // запись температуры
		    for (i=0;i<nodes; i++) {
			    //fprintf(fp, "%e ", potent[TEMP][i]);
			    fprintf(fp, "%e ", potent[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");

		}
		/*
		// запись горизонтальной скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
        // запись вертикальной скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

		// запись модуля скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        */

		if (btriangle) {
			// каждый четырёхугольник состоит из двух треугольников.
			for (i=0;i<nelmts; i++) {
				fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[3][i]);
				fprintf(fp, "%d %d %d\n", nvtx[1][i], nvtx[3][i], nvtx[2][i]);
			}
		}
		else {
		    // запись информации о разностной сетке
		    for (i=0;i<nelmts; i++) {
			    if (nve==3) fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i]);
			    if (nve==4) fprintf(fp, "%d %d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i]);
		    }
		}

		fclose(fp); // закрытие файла
        //WinExec("C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2008\\bin\\tec360.exe temp_xy.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2009\\bin\\tec360.exe temp_xy.PLT", SW_NORMAL);
		/*
		STARTUPINFO cif;
        ZeroMemory(&cif,sizeof(STARTUPINFO));
		cif.cb = sizeof(STURTUPINFO);
        PROCESS_INFORMATION pi;
		ZeroMemory(&pi,sizeof(PROCESS_INFORMATION));
 
        
        CreateProcess(NULL, L"C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2009\\bin\\tec360.exe temp_xy.PLT",
        NULL,
        NULL,NULL,FALSE,NULL,NULL,NULL,&cif,&pi);
		*/


		if (btriangle) {
			// граничные условия
			FILE *fp1;
	        errno_t err1;
            if ((err1 = fopen_s( &fp1, "gru.txt", "r")) != 0) {
		         printf("Open gru.txt File Error\n");
	        }
	        else {

				FILE *fp2;
	            errno_t err2;
	            // создание файла для записи.
	            if ((err2 = fopen_s( &fp2, "tigel_boundary_tec.dat", "w")) != 0) {
		            printf("Create File Error\n");
	            }
	            else {

				    int inumbound=0;
				    float fiso, fstart, fend, fval;
					int iorient=0; // 0-x, 1-y.

				    fscanf_s(fp1, "%d", &inumbound);

					fprintf(fp2,"BOUND = %d\n",inumbound);
					fprintf(fp2,"TITLE = \"mesh\"\n");
					fprintf(fp2,"VARIABLES = x, y\n");

					for (int j=0; j<inumbound; j++) {
					     fscanf_s(fp1, "%d", &iorient);
						 fscanf_s(fp1, "%f", &fiso);
						 fscanf_s(fp1, "%f", &fstart);
						 fscanf_s(fp1, "%f", &fend);
						 fscanf_s(fp1, "%f", &fval);

						 int inod=0;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 inod++;
									 }
								 }
							 }

		                 }

			             fprintf(fp2,"ZONE T=\"edge%d\", N=%d, E=%d, ET=LINESEG, F=FEBLOCK, VAL=%1.4f\n",j,inod,inod-1,fval);
						 fprintf(fp2,"\n");
						 inod=0;
						 // печать x
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }

		                 }
						 fprintf(fp2, "\n");
						  fprintf(fp2, "\n");

						  inod=0;
						 // печать y
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }

		                 }
						 fprintf(fp2, "\n");
						  fprintf(fp2, "\n");

						 // пары должны иметь глобальную нумерацию начиная с единицы.
						 int iprev=-1;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 if (iprev>-1) {
											 fprintf(fp2,"%d %d \n",iprev,i+1);
										 }
										 iprev=i+1;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 if (iprev>-1) {
											 fprintf(fp2,"%d %d \n",iprev,i+1);
										 }
										 iprev=i+1;
									 }
								 }
							 }

		                 }


					}

					fclose(fp2);

	            }

				fclose(fp1);
	        }
		}

	}
} // exporttecplotxy360

// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// универсально подходит и для треугольных и для квадратных ячеек.
void exporttecplotxy360extrude3D(int nve, int nelmts, int nodes, int** nvtx, Real* x, Real* y, Real* potent)
{
	FILE *fp;
	errno_t err;

	float hz=0.2; // шаг в третьем измерении.
	int nz=10; // количество слоёв в 3D (не менее 2ух).

	// создание файла для записи.
	if ((err = fopen_s( &fp, "tigel_interior_tec.dat"/*"temp_xy.PLT"*/, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// запись заголовка
		fprintf(fp, "TITLE = \"temperature\"\n");

		// запись имён переменных
		//fprintf(fp, "VARIABLES = x, y, \"temperature\" , Vx, Vy, Mag \n");
		if (btriangle) {
			fprintf(fp, "VARIABLES = x, y, z \n"); // без температуры (зачем лишние данные это же сеточный генератор).
		}
		else {
		   fprintf(fp, "VARIABLES = x, y, z, \"temperature\" \n");
		}

		// запись информации о зонах
		if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nz*nodes, nelmts*(nz-1));
        if (nve==4) {
			if (btriangle) {
				// каждый четырёхугольник это два треугольника, поэтому 
				// треугольных элементов в два раза больше.
				// Дело в том что поддерживается только BRICK и TETRAHEDRON
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=PRIZM, F=FEBLOCK\n\n", nz*nodes, 2*nelmts*(nz-1));
			}
			else {
				// истинные четырёхугольные элементы.
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", nz*nodes, nelmts*(nz-1));
			}
		}

		int i=0; // счётчики 
		int j=0; // цикла for

		// запись x
		for (int k=0; k<nz; k++) {
	    for (i=0; i<nodes; i++) {	
				fprintf(fp, "%e ", x[i]);
				if (i%10==0) fprintf(fp, "\n");
		}
		}
			
		fprintf(fp, "\n");
          
		// запись y
		for (int k=0; k<nz; k++) {
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", y[i]);
            if (i%10==0) fprintf(fp, "\n");
		}
		}
			
        fprintf(fp, "\n");

		// запись z
		for (int k=0; k<nz; k++) {
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", 0.0+k*hz);
            if (i%10==0) fprintf(fp, "\n");
		}
		}

		if (!btriangle) {
		    // запись температуры
			for (int k=0; k<nz; k++) {
		    for (i=0;i<nodes; i++) {
			    //fprintf(fp, "%e ", potent[TEMP][i]);
			    fprintf(fp, "%e ", potent[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }
			}

		    fprintf(fp, "\n");

		}
		/*
		// запись горизонтальной скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
        // запись вертикальной скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

		// запись модуля скорости
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        */

		if (btriangle) {
			// каждый четырёхугольник состоит из двух треугольников.
			for (int k=0; k<nz-1; k++) {
			for (i=0;i<nelmts; i++) {
				fprintf(fp, "%d %d %d %d %d %d\n", nvtx[0][i]+k*nodes, nvtx[1][i]+k*nodes, nvtx[3][i]+k*nodes,nvtx[0][i]+(k+1)*nodes, nvtx[1][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes);
				fprintf(fp, "%d %d %d %d %d %d\n", nvtx[1][i]+k*nodes, nvtx[3][i]+k*nodes, nvtx[2][i]+k*nodes,nvtx[1][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes, nvtx[2][i]+(k+1)*nodes);
			}
			}
		}
		else {
		    // запись информации о разностной сетке
			for (int k=0; k<nz-1; k++) {
		    for (i=0;i<nelmts; i++) {
			    if (nve==3) fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i]);
			    if (nve==4) {
					fprintf(fp, "%d %d %d %d %d %d %d %d\n", nvtx[0][i]+k*nodes, nvtx[1][i]+k*nodes, nvtx[2][i]+k*nodes, nvtx[3][i]+k*nodes, nvtx[0][i]+(k+1)*nodes, nvtx[1][i]+(k+1)*nodes, nvtx[2][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes);
				}
		    }
			}
		}

		fclose(fp); // закрытие файла
        //WinExec("C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2008\\bin\\tec360.exe temp_xy.PLT",SW_NORMAL);
		  //WinExec("C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2009\\bin\\tec360.exe temp_xy.PLT", SW_NORMAL);
		/*
		STARTUPINFO cif;
        ZeroMemory(&cif,sizeof(STARTUPINFO));
		cif.cb = sizeof(STURTUPINFO);
        PROCESS_INFORMATION pi;
		ZeroMemory(&pi,sizeof(PROCESS_INFORMATION));
 
        
        CreateProcess(NULL, L"C:\\Program Files\ (x86)\\Tecplot\\Tec360\ 2009\\bin\\tec360.exe temp_xy.PLT",
        NULL,
        NULL,NULL,FALSE,NULL,NULL,NULL,&cif,&pi);
		*/

		
		/*if (btriangle) {
			// граничные условия
			FILE *fp1;
	        errno_t err1;
            if ((err1 = fopen_s( &fp1, "gru.txt", "r")) != 0) {
		         printf("Open gru.txt File Error\n");
	        }
	        else {

				FILE *fp2;
	            errno_t err2;
	            // создание файла для записи.
	            if ((err2 = fopen_s( &fp2, "tigel_boundary_tec.dat", "w")) != 0) {
		            printf("Create File Error\n");
	            }
	            else {

				    int inumbound=0;
				    float fiso, fstart, fend, fval;
					int iorient=0; // 0-x, 1-y.

				    fscanf_s(fp1, "%d", &inumbound);

					fprintf(fp2,"BOUND = %d\n",inumbound);
					fprintf(fp2,"TITLE = \"mesh\"\n");
					fprintf(fp2,"VARIABLES = x, y\n");

					for (int j=0; j<inumbound; j++) {
					     fscanf_s(fp1, "%d", &iorient);
						 fscanf_s(fp1, "%f", &fiso);
						 fscanf_s(fp1, "%f", &fstart);
						 fscanf_s(fp1, "%f", &fend);
						 fscanf_s(fp1, "%f", &fval);

						 int inod=0;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 inod++;
									 }
								 }
							 }

		                 }

			             fprintf(fp2,"ZONE T=\"edge%d\", N=%d, E=%d, ET=LINESEG, F=FEBLOCK, VAL=%1.4f\n",j,inod,inod-1,fval);
						 fprintf(fp2,"\n");
						 inod=0;
						 // печать x
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }

		                 }
						 fprintf(fp2, "\n");
						  fprintf(fp2, "\n");

						  inod=0;
						 // печать y
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }

		                 }
						 fprintf(fp2, "\n");
						  fprintf(fp2, "\n");

						 // пары должны иметь глобальную нумерацию начиная с единицы.
						 int iprev=-1;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // горизонтальная линия.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 if (iprev>-1) {
											 fprintf(fp2,"%d %d \n",iprev,i+1);
										 }
										 iprev=i+1;
									 }
								 }
							 }
							 if (iorient==1) {
								 // вертикальная линия.
								 if (fabs(x[i]-fiso)<1.0e-20) {
									 if ((y[i]>fstart)&&(y[i]<fend)) {
										 if (iprev>-1) {
											 fprintf(fp2,"%d %d \n",iprev,i+1);
										 }
										 iprev=i+1;
									 }
								 }
							 }

		                 }


					}

					fclose(fp2);

	            }

				fclose(fp1);
	        }
		}*/

	}
} // exporttecplotxy360extrude3D


// конец


