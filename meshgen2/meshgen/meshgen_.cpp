// meshgen.cpp: определяет точку входа для консольного приложения.
//  Это консольный сеточный генератор не требующий вмешательства пользователя.
// Основан на алгебраических методах.
// Dynamic Algebraic Mesh Generator v0.0.1 30 ноября 2012.

#include "stdafx.h"
#include <stdio.h>
#include <tchar.h>
#include <conio.h>
#include <math.h>

#include "AliceDataDubl.cpp"


const double eps=1e-20; // для отделения вещественного нуля.

#define NONE -1 // прямая линия
#define MYPARABOLA 0 // парабола
#define MYSIN 1 // синусоида
#define MYCIRCLE 2 // окружность

#define CURVETOP 0 // NONE -1 просто прямоугольник.
#define CURVEBOTTOM 1
#define CURVEDOUBLE 2

typedef enum States { Normal, Slash, Comment,EndComment, EndFile,StrState,SimState,OneComment /* перечисление всех
                                                  состояний автомата */ } States;

typedef struct TBLOCK {

	int priority; // приоритет блока.

	double xS, yS, xE, yE; // размеры блока
	int inhollow; // если 0 то hollow блок.
	double rho, cp, lambda; // свойства материала.
} BLOCK;

// Универсальный блок какие либо две противоположные границы
// этого элемента задаётся функцией.
typedef struct TBLOCKU {

	int priority; // приоритет блока.

	int inx, iny; // число узлов сетки по горизонтали и вертикали.

	// yE граница верхней зоны, этот элемент должен
	// быть вписан в прямоугольник чтобы всё было согласованно с другими блоками.

	int icurvetop, icurvebottom; // выбор кривой сверху и снизу.

	double xS, yS, xE, yE; // размеры блока

	int itoporbottom; // где задаются кривые сверху или снизу.


	// t - top, b - bottom.
	double at,bt,ct; // формула для параболы : a*x^2+b*x+c (icurve==MYPARABOLA)
	double dyt, epsat, omegat, dxt; // формула для синуса : dy+epsa*sin(2*MPI*omega*(x-dx)) (icurve==MYSIN).
	double ab,bb,cb; // формула для параболы : a*x^2+b*x+c (icurve==MYPARABOLA)
	double dyb, epsab, omegab, dxb; // формула для синуса : dy+epsa*sin(2*MPI*omega*(x-dx)) (icurve==MYSIN).

	// в случае icurvetop==MYCIRCLE верхняя граница задаётся в виде дуги окружности.
	double yCt, rt, angleSt, angleEt; // y=yCt+rt*sin(angle).
	
	int inhollow; // если 0 то hollow блок.
	double rho, cp, lambda; // свойства материала.
} BLOCKU;

// универсальный сектор
typedef struct TSECTORU {
	int priority; // приоритет блока.

	int inr, intheta; // число узлов сетки по горизонтали и вертикали.
	int ip; // идентификатор перехода с круга на квадрат и наоборот.
	// что обозначают значения ip перечислено в файле obj.txt.

	double xC, yC; // центр окружности.
	double rmin, rmax; // диапазон изменения радиуса.
	double thetamin, thetamax; // диапазон изменения угла.

	int inhollow; // если 0 то hollow блок.
	double rho, cp, lambda; // свойства материала.

} SECTORU;

typedef struct TBLOCKNEW {
	int priority; // приоритет блока.

	int n1, n2; // количество разбиений по сторонам (1,3) и (2,4)

	double x1, y1, x2, y2, x3, y3, x4, y4;

	int inhollow; // если 0 то hollow блок.
	double rho, cp, lambda; // свойства материала.
} BLOCKNEW;

BLOCK* b; // прямоугольные блоки.
BLOCKU* bt;
SECTORU* bsec; // сектор.
BLOCKNEW* bnew; // произвольный выпуклый четырёхугольник.

// минимальное количество клеток по 
// каждому из направлений на минимальном
// элементе разбиения
int min_elem_in_x_element=4; // обязательно не меньше 4!
int min_elem_in_y_element=4; // обязательно не меньше 4!
// ручная адаптация сетки
int adapt_x=0;
int adapt_y=0;

// реализация функции SetLength из Delphi 
// для изменения размеров динамического массива ra. 
// Возвращает модифицированный массив ra необходимого 
//  размера. 
// Во время работы функции утечки памяти не происходит.
void SetLength(double* &ra, int isizeold, int isize) 
{

	// isize - новая длина динамического массива.
	double *temp;
	temp = new double[isize];
    int i;
	for (i=0; i<isize; i++) temp[i]=0.0; // инициализация
    /*
	int isizeold;
	if (ra != NULL) {
       isizeold = sizeof(ra)/sizeof(ra[0]); // длина старого массива
       printf("%d  ",isizeold); // не хочет правильно определять размер массива 
	} 
	  else
	{
       isizeold=0;
	}
	*/
	if (isize < isizeold) isizeold=isize;
	for (i=0; i<isizeold; i++) temp[i]=ra[i];
	
	
	if (isizeold!=0) {
		delete  ra; // уничтожение объекта
	}
	ra = new double[isize]; // выделение памяти
	for (i=0; i<isize; i++) ra[i]=temp[i]; // копирование
	delete temp; // освобождение памяти
	
} // SetLength

// добавляет несуществующую границу к массиву
void addboundary(double* &rb, int &in, double g) {
	// rb - модифицируемый массив границ,
	// in - номер последней границы в массиве, нумерация начинается с нуля.
	// g - граница, кандидат на добавление.
	int i=0;
    
	bool bfind=false;
	for (i=0; i<=in; i++) if (fabs(rb[i]-g)<eps) bfind=true;
	if (!bfind) {
        SetLength(rb, in+1, in+2);
		in++;
		rb[in]=g;
	}
} // addboundary


// Сортировка по возрастанию прямым обменом -
// пузырьковая улучшенная.
// Основной принцип: если элементы уже упорядочены, 
// сортировка прекращается.
// in - предполагается достаточно малым меньше 500,
// именно поэтому применима пузырьковая сортировка.
void BubbleEnhSort(double* &rb, int in) {
	int i,j,k;
	double x;
	bool swapped; // факт обмена

	for (i=1; i<=in; i++) {
		k=0; // количество обменов данного прохода
		swapped=false; // обменов не было
		// i пузырьков уже всплыло, их просматривать ненужно
		for (j=in; j>=i; j--) {
			if (rb[j-1]>rb[j]) {
				// SWAP
				x=rb[j-1];
				rb[j-1]=rb[j];
				rb[j]=x;
				k++;
				swapped=true;
			}
		}
		if (!swapped) break; // выход из цикла
	}
} // BubbleEnhSort

// Целочисленная сортировка по возрастанию прямым обменом -
// пузырьковая улучшенная.
// Основной принцип: если элементы уже упорядочены, 
// сортировка прекращается.
// in - предполагается достаточно малым меньше 200,
// именно поэтому применима пузырьковая сортировка.
void BubbleEnhSorti(int* &rb, int in) {
	int i,j,k;
	int x;
	bool swapped; // факт обмена

	for (i=1; i<=in; i++) {
		k=0; // количество обменов данного прохода
		swapped=false; // обменов не было
		// i пузырьков уже всплыло, их просматривать ненужно
		for (j=in; j>=i; j--) {
			if (rb[j-1]>rb[j]) {
				// SWAP
				x=rb[j-1];
				rb[j-1]=rb[j];
				rb[j]=x;
				k++;
				swapped=true;
			}
		}
		if (!swapped) break; // выход из цикла
	}
} // BubbleEnhSorti

double fmin(double ra, double rb) {
	double r=ra;
	if (rb<ra) r=rb;
	return rb;
} // fmin

// возвращает наибольшее из двух чисел
double fmax(double fA, double fB) {
	double r=fA;
	if (fB>r) r=fB;
	return r;
} // fmax

/* генерирует простейшую расчётную сетку hex cartesian - 
 * структурированная прямоугольная сетка в трёхмерной расчётной 
 * области в виде кирпича.
 * В функцию передаются значения общего количества делений по каждой из осей inx, iny. 
 * Число блоков прямоугольной формы и размеры этих блоков.
 * В результате получаются массивы xpos и ypos содержащие позиции сеточных узлов.
*/
void simplemeshgen(double* &xpos, double* &ypos, int &inx, int &iny, 
				   int lb, BLOCK* b)
{

	bool bgeom=true; // если true, то используется неравномерная сетка по закону геометрической прогрессии.
	double q=1.2; // 1.05 - очень слабо неравномерная, 1.1, 1.2 - умеренно неравномерная, 1.25, 1.5, 2 - сильно неравномерная. 
	// для низкорейнольдсовых моделей турбулентности знаменатель геометрической прогрессии должен быть не больше 1.3.

	int i;
	// По Оси Ox
	double *rxboundary; // массив обязательных границ
    int inumboundaryx=1; 
	rxboundary = new double [inumboundaryx+1]; // число границ по оси x

	// добавить границы кабинета
	rxboundary[0]=b[0].xS; // начало области
	rxboundary[inumboundaryx]=b[0].xE; // конец области
	
	// блоки
	for (i=1; i<lb; i++) {
		addboundary(rxboundary,inumboundaryx,b[i].xS);
        addboundary(rxboundary,inumboundaryx,b[i].xE);
	}

    // упорядочивание по возрастанию
	BubbleEnhSort(rxboundary, inumboundaryx);

    int *ixintervalcount; // число интервалов
	ixintervalcount = new int [inumboundaryx]; // на один меньше чем число границ.
	double alphascal;
	int inowintervalcount;
	for (i=0; i<(inumboundaryx); i++) {
         alphascal=(rxboundary[i+1]-rxboundary[i])/(rxboundary[inumboundaryx]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*inx);
		 if (inowintervalcount < min_elem_in_x_element) inowintervalcount=min_elem_in_x_element;
         ixintervalcount[i]=inowintervalcount;
	}
    /* // debug
	for (i=0; i<(inumboundaryx); i++) {
         printf("%d  ",ixintervalcount[i]);
	} //*/
    //getchar();

	// заполнение масива положения узлов.
	int iposmark = 1;
	double dx;
	int k;
	int ixoldsize=0;
	SetLength(xpos,ixoldsize,1);
    ixoldsize=1;
	for (i=0; i<(inumboundaryx); i++) 
	{
		if ((ixintervalcount[i]-2)%2==0) ixintervalcount[i]++; // чтобы число узлов было всегда чётное.
		int n2=(int)((ixintervalcount[i]-1)/2); // округление в меньшую сторону
		double qn2=q;
		for (int i1=1; i1<n2; i1++) qn2*=q; // возведение в степень.
		double b1=(rxboundary[i+1]-rxboundary[i])*(q-1.0)/(2.0*(qn2-1.0));
		//printf("length=%e\n",(rxboundary[i+1]-rxboundary[i]));
		//getchar(); // OK

        dx=(rxboundary[i+1]-rxboundary[i])/(ixintervalcount[i]-1);
		SetLength(xpos,ixoldsize,(iposmark+ixintervalcount[i]));
        ixoldsize=(iposmark+ixintervalcount[i]);
        //printf("%d  ",ixoldsize);// debug
        for (k=iposmark; k<=iposmark+ixintervalcount[i]-2; k++)
		{
			if (!bgeom) {
				// равномерная сетка
				xpos[k]=rxboundary[i]+(k-iposmark)*dx;
			}
			else
			{
				// неравномерная сетка
				int ic1=k-iposmark;
				double madd=b1;
				if (ic1<=n2) {
					// сначала левая сторона.
					for (int i1=1; i1<ic1; i1++) madd*=q; 
				} else 
				{
					// потом правая сторона.
					for (int i1=2*n2-ic1; i1>0; i1--) madd*=q;
				}
				if (k==iposmark) xpos[k]=rxboundary[i];
				else xpos[k]=xpos[k-1]+madd;
				
			}
		}
		iposmark=iposmark+ixintervalcount[i]-1;
	}
	SetLength(xpos,ixoldsize,iposmark+1);
	xpos[iposmark]=rxboundary[inumboundaryx];
	inx=iposmark;
    for (i=0; i<inx; i++) xpos[i]=xpos[i+1]; // сдвиг влево на 1
    SetLength(xpos,inx+1,inx); 
	inx--; // индексация начинается с нуля и заканчивается значением inx
	
	//for (i=0; i<adapt_x; i++) simplecorrect_meshgen_x(xpos, inx, lb, ls, lw, b, s, w);

	/* // debug
	for (i=0; i<=inx; i++) {
        printf("%f  ",xpos[i]);   
	}
	getchar(); //*/ 


    // По оси Oy
    double *ryboundary; // массив обязательных границ
    int inumboundaryy=1;
	ryboundary = new double [inumboundaryy+1]; // число границ по оси y
	ryboundary[0]= b[0].yS; // начало области 
	ryboundary[inumboundaryy]= b[0].yE; // конец области 
    
    // блоки
	for (i=1; i<lb; i++) {
		addboundary(ryboundary,inumboundaryy,b[i].yS);
        addboundary(ryboundary,inumboundaryy,b[i].yE);
	}


    // упорядочивание по возрастанию
	BubbleEnhSort(ryboundary, inumboundaryy);


    int *iyintervalcount; // число интервалов
    iyintervalcount = new int [inumboundaryy]; // на один меньше чем число границ.
    for (i=0; i<(inumboundaryy); i++) {
         alphascal=(ryboundary[i+1]-ryboundary[i])/(ryboundary[inumboundaryy]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*iny);
		 if (inowintervalcount < min_elem_in_y_element) inowintervalcount=min_elem_in_y_element; 
         iyintervalcount[i]=inowintervalcount;
	}
    // заполнение массива узлов
    iposmark = 1;
	double dy;
    int iyoldsize=0;
    SetLength(ypos,iyoldsize,1);
    iyoldsize=1;
	for (i=0; i<inumboundaryy; i++) 
	{

		// чтобы число узлов було всегда чётное
		if ((iyintervalcount[i]-2)%2==0) iyintervalcount[i]++;
		int n2=(int)((iyintervalcount[i]-1)/2); // округление в меньшую сторону
		double qn2=q;
		for (int i1=1; i1<n2; i1++) qn2*=q; // возведение в степень.
		double b1=(ryboundary[i+1]-ryboundary[i])*(q-1.0)/(2.0*(qn2-1.0));

        dy=(ryboundary[i+1]-ryboundary[i])/(iyintervalcount[i]-1);
		SetLength(ypos,iyoldsize,iposmark+iyintervalcount[i]);
        iyoldsize=iposmark+iyintervalcount[i];
        for (k=iposmark; k<=iposmark+iyintervalcount[i]-2; k++) 
		{
			if (!bgeom) {
				// равномерная сетка
			    ypos[k]=ryboundary[i]+(k-iposmark)*dy;
			}
			else 
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				int ic1=k-iposmark;
				double madd=b1;
				if (ic1<=n2) {
					for (int i1=1; i1<ic1; i1++) madd*=q;
				} else 
				{
					for (int i1=2*n2-ic1; i1>0; i1--) madd*=q;
				}
				if (k==iposmark) ypos[k]=ryboundary[i];
				else ypos[k]=ypos[k-1]+madd;
			}
		}
		iposmark=iposmark+iyintervalcount[i]-1;
	}
	SetLength(ypos,iyoldsize,iposmark+1);
	ypos[iposmark]=ryboundary[inumboundaryy];
	iny=iposmark;
    for (i=0; i<iny; i++) ypos[i]=ypos[i+1]; // сдвиг влево на 1
    SetLength(ypos,iny+1,iny); 
	iny--; // индексация начинается с нуля и заканчивается значением iny
    
    //for (i=0; i<adapt_y; i++) simplecorrect_meshgen_y(ypos, iny, lb, ls, lw, b, s, w);


	/* // debug
	for (i=0; i<=iny; i++) {
        printf("%f  ",ypos[i]);   
	}//*/


    // Освобождение оперативной памяти.
	delete rxboundary;
	delete ryboundary;
	
	delete ixintervalcount;
	delete iyintervalcount;
    
} // simplemeshgen

// удаляет коментарии  из входного файла.
int replace() {
	 FILE  *fi,  *fo;         /* входной и выходной файл */
     States State = Normal;   /* текущее состояние */
     int/*space*/c;                   /* считанный символ (только один текущий!) */
     /* все, никакие переменные больше не нужны, этих вполне достаточно */
     errno_t err;
	
  
  if ((err = fopen_s(&fi,"obj.txt", "rb"))!=0) {
    fprintf(stderr, "Input file obj.txt open error.\n");
    return 1;
  }
  if ((err = fopen_s(&fo,"no_coment.txt", "wb"))!=0) {
    fclose(fi);
    fprintf(stderr, "Output file no_coment.txt open error.\n");
    return 2;
  }
 
  while ((c=fgetc(fi)) != EOF)  /* считываем символ и проверяем, не конец ли файла? */
   {
    switch (State)   /* если нет, то обрабатываем в зависимости от текущего состояния */
    {
      case Normal:
 
                        if (c == '/')         /* если встретили слэш, */
                                 State = Slash;     /* то перешли в соответствующее состояние */
                        else {
                                 fputc(c,fo);
                                 if (c == '"') State=StrState;
                                 if (c == '\'') State=SimState;
                              }
 
                        break;
 
        case Slash:          /* если предыдущим был слэш, то ... */
            if (c == '*')
                State=Comment;
            else if (c == '/') 
                State=OneComment;
                else              
                 
            {
                State=Normal;
                fputc('/',fo);
                fputc(c,fo);
            }
             break;
 
        case Comment:
            if ( c=='*')
                State=EndComment;
            //else State=Comment;
            break;
 
        case EndComment:
            State=Comment;
            if (c=='/')
                {State=Normal; fputc(' ',fo);}
            if (c=='*')
                State=EndComment;
            break;
 
        case StrState:
            fputc(c,fo);
            if (c=='"') State= Normal; else
            if (c == '\\')
            {
            if ((c=fgetc(fi)) == EOF) break;
            fputc(c,fo);
            }
        break;
 
        case SimState:
            fputc(c,fo);
            if (c=='\'') State= Normal; else
            if (c == '\\')
            {
            if ((c=fgetc(fi)) == EOF) break;
            fputc(c,fo);
            }
        break;
 
        case OneComment:
                 if ( (c == '\r') || (c=='\n'))
                 {
                 State=Normal;
                 fputc(c,fo);}
        break;
        
        }
   }      if (State==Slash)   fputc(c,fo);
 
  fclose(fi);    /* не забывайте закрывать файлы! */
  fclose(fo);    /* особенно выходной, открытый для записи */
  return 0;       /* нормальное завершение работы функции*/
} // replace 

// считывает объекты из которых состоит сетка.
// lb - число блоков.
int readobj(int &inx, int &iny, int &lb, int &lbt, int &lbsec, int &lbnew) {

	FILE *fp;
    errno_t err;
    if ((err = fopen_s(&fp,"no_coment.txt", "r"))!=0) {
       fprintf(stderr, "Input file no_coment.txt open error.\n");
	   return 1;
    }

    float fin=0.0;
	int din=0;
	int imaxblock=0;

	fscanf_s(fp, "%d", &din);
	inx=din;
	fscanf_s(fp, "%d", &din);
	iny=din;

	fscanf_s(fp, "%d", &din);
	imaxblock=din;
	lb=imaxblock;

	fscanf_s(fp, "%d", &din);
	lbt=din; // число блоков с верхней или нижней границей заданной в виде параболы или синусоиды.

	fscanf_s(fp, "%d", &din);
	lbsec=din; // число секторов.

	// Количество новых произвольных выпуклых многоугольников.
    fscanf_s(fp, "%d", &din);
	lbnew=din; 

	b=new BLOCK[imaxblock]; // выделение памяти под прямоугольные блоки.
	bt=new BLOCKU[lbt]; // выделение памяти под блоки у которых верхняя граница задана в виде параболы.
	bsec=new SECTORU[lbsec]; // выделение памяти под сектора. 
	bnew=new BLOCKNEW[lbnew]; // произвольный выпуклый четырёхугольник.

	int ib=0, ibt=0, ibsec=0, ibnew=0;

	for (int i=0; i<lb+lbt+lbsec+lbnew; i++) {

		int id;
		fscanf_s(fp, "%d", &din);
	    id=din;


		if (id==0) {
			// прямоугольный блок.

            fscanf_s(fp, "%d", &din);
			b[ib].inhollow=din; // если 0 то hollow block.

			fscanf_s(fp, "%d", &din);
			b[ib].priority=din; // приоритет блока.

			// размеры блока.
	        fscanf_s(fp, "%f", &fin);
			b[ib].xS=fin;
            fscanf_s(fp, "%f", &fin);
			b[ib].yS=fin;
			 fscanf_s(fp, "%f", &fin);
			b[ib].xE=fin;
			 fscanf_s(fp, "%f", &fin);
			b[ib].yE=fin;

			// свойства материала.
			fscanf_s(fp, "%f", &fin);
			b[ib].rho=fin;
			fscanf_s(fp, "%f", &fin);
			b[ib].cp=fin;
			fscanf_s(fp, "%f", &fin);
			b[ib].lambda=fin;

			ib++;
		}

		if (id==1) {
			// прямоугольный блок.

            fscanf_s(fp, "%d", &din);
			bt[ibt].inhollow=din; // если 0 то hollow block.

			fscanf_s(fp, "%d", &din);
			bt[ibt].priority=din; // приоритет блока.

			fscanf_s(fp, "%d", &din);
			bt[ibt].inx=din;
			fscanf_s(fp, "%d", &din);
			bt[ibt].iny=din;

			// размеры блока.
	        fscanf_s(fp, "%f", &fin);
			bt[ibt].xS=fin;
            fscanf_s(fp, "%f", &fin);
			bt[ibt].yS=fin;
			 fscanf_s(fp, "%f", &fin);
			bt[ibt].xE=fin;
			 fscanf_s(fp, "%f", &fin);
			bt[ibt].yE=fin;

			fscanf_s(fp, "%d", &din);
			bt[ibt].icurvetop=din;

			fscanf_s(fp, "%d", &din);
			bt[ibt].icurvebottom=din;

			fscanf_s(fp, "%d", &din);
			bt[ibt].itoporbottom=din;

			if (bt[ibt].icurvetop==NONE) {
				bt[ibt].at=0.0;
				bt[ibt].bt=0.0;
				bt[ibt].ct=0.0;
				bt[ibt].dyt=0.0;
			    bt[ibt].epsat=0.0;
			    bt[ibt].omegat=0.0;
			    bt[ibt].dxt=0.0;

				bt[ibt].yCt=0.0;
				bt[ibt].rt=0.0;
                bt[ibt].angleSt=0.0;
				bt[ibt].angleEt=0.0;
			}
			else if (bt[ibt].icurvetop==MYPARABOLA) {
			   // парабола.
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].at=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].bt=fin;
			   fscanf_s(fp, "%f", &fin);

			   bt[ibt].ct=fin;
			   bt[ibt].dyt=0.0;
			   bt[ibt].epsat=0.0;
			   bt[ibt].omegat=0.0;
			   bt[ibt].dxt=0.0;

			    bt[ibt].yCt=0.0;
				bt[ibt].rt=0.0;
                bt[ibt].angleSt=0.0;
				bt[ibt].angleEt=0.0;

			}
			else if (bt[ibt].icurvetop==MYSIN) {
				 fscanf_s(fp, "%f", &fin);
			     bt[ibt].dyt=fin;
				 fscanf_s(fp, "%f", &fin);
			     bt[ibt].epsat=fin;
			     fscanf_s(fp, "%f", &fin);
			     bt[ibt].omegat=fin;
			     fscanf_s(fp, "%f", &fin);

			     bt[ibt].dxt=fin;
				 bt[ibt].at=0.0;
                 bt[ibt].bt=0.0;
				 bt[ibt].ct=0.0;

				 bt[ibt].yCt=0.0;
				 bt[ibt].rt=0.0;
                 bt[ibt].angleSt=0.0;
				 bt[ibt].angleEt=0.0;

			}
			else if (bt[ibt].icurvetop==MYCIRCLE) {
				 fscanf_s(fp, "%f", &fin);
				 bt[ibt].yCt=fin;
				 fscanf_s(fp, "%f", &fin);
				 bt[ibt].rt=fin;
			     fscanf_s(fp, "%f", &fin);
				 bt[ibt].angleSt=fin;
			     fscanf_s(fp, "%f", &fin);
				 bt[ibt].angleEt=fin;
				                

				 bt[ibt].dyt=0.0;
				 bt[ibt].epsat=0.0;
                 bt[ibt].omegat=0.0;
				 bt[ibt].dxt=0.0;

				 bt[ibt].at=0.0;
                 bt[ibt].bt=0.0;
				 bt[ibt].ct=0.0;
			}
		    else {
				// парабола. (по умолчанию)
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].at=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].bt=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].ct=fin;
			   bt[ibt].dyt=0.0;
			   bt[ibt].epsat=0.0;
			   bt[ibt].omegat=0.0;
			   bt[ibt].dxt=0.0;

			   bt[ibt].yCt=0.0;
			   bt[ibt].rt=0.0;
               bt[ibt].angleSt=0.0;
			   bt[ibt].angleEt=0.0;

			   bt[ibt].icurvetop=MYPARABOLA;
			}

			if (bt[ibt].icurvebottom==NONE) {
				bt[ibt].ab=0.0;
				bt[ibt].bb=0.0;
				bt[ibt].cb=0.0;
				bt[ibt].dyb=0.0;
			    bt[ibt].epsab=0.0;
			    bt[ibt].omegab=0.0;
			    bt[ibt].dxb=0.0;
			}
			else if (bt[ibt].icurvebottom==MYPARABOLA) {
			   // парабола.
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].ab=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].bb=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].cb=fin;
			   bt[ibt].dyb=0.0;
			   bt[ibt].epsab=0.0;
			   bt[ibt].omegab=0.0;
			   bt[ibt].dxb=0.0;
			}
			else if (bt[ibt].icurvebottom==MYSIN) {
				 fscanf_s(fp, "%f", &fin);
			     bt[ibt].dyb=fin;
				 fscanf_s(fp, "%f", &fin);
			     bt[ibt].epsab=fin;
			     fscanf_s(fp, "%f", &fin);
			     bt[ibt].omegab=fin;
			     fscanf_s(fp, "%f", &fin);
			     bt[ibt].dxb=fin;
				 bt[ibt].ab=0.0;
                 bt[ibt].bb=0.0;
				 bt[ibt].cb=0.0;
			}
			else {
				// парабола. (по умолчанию)
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].ab=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].bb=fin;
			   fscanf_s(fp, "%f", &fin);
			   bt[ibt].cb=fin;
			   bt[ibt].dyb=0.0;
			   bt[ibt].epsab=0.0;
			   bt[ibt].omegab=0.0;
			   bt[ibt].dxb=0.0;
			   bt[ibt].icurvebottom=MYPARABOLA;
			}

			// свойства материала.
			fscanf_s(fp, "%f", &fin);
			bt[ibt].rho=fin;
			fscanf_s(fp, "%f", &fin);
			bt[ibt].cp=fin;
			fscanf_s(fp, "%f", &fin);
			bt[ibt].lambda=fin;

			ibt++;
		}

		if (id==2) {
			// сектор.

            fscanf_s(fp, "%d", &din);
			bsec[ibsec].inhollow=din; // если 0 то hollow block.

			fscanf_s(fp, "%d", &din);
			bsec[ibsec].priority=din; // приоритет блока.

			fscanf_s(fp, "%d", &din);
			bsec[ibsec].inr=din;
			fscanf_s(fp, "%d", &din);
			bsec[ibsec].intheta=din;

			// размеры блока.
	        fscanf_s(fp, "%f", &fin);
			bsec[ibsec].xC=fin;
            fscanf_s(fp, "%f", &fin);
			bsec[ibsec].yC=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].rmin=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].thetamin=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].rmax=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].thetamax=fin;

			fscanf_s(fp, "%d", &din);
			bsec[ibsec].ip=din; // идентификатор перхода. 

			// свойства материала.
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].rho=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].cp=fin;
			fscanf_s(fp, "%f", &fin);
			bsec[ibsec].lambda=fin;

			ibsec++;
		}

		if (id==3) {
			// произвольный выпуклый четырёхугольник.
			fscanf_s(fp, "%d", &din);
			bnew[ibnew].inhollow=din; // если 0 то hollow block.

			fscanf_s(fp, "%d", &din);
			bnew[ibnew].priority=din; // приоритет блока.

			fscanf_s(fp, "%d", &din);
			bnew[ibnew].n1=din;
			fscanf_s(fp, "%d", &din);
			bnew[ibnew].n2=din;

			// размеры блока.
	        fscanf_s(fp, "%f", &fin);
			bnew[ibnew].x1=fin;
            fscanf_s(fp, "%f", &fin);
			bnew[ibnew].y1=fin;
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].x2=fin;
            fscanf_s(fp, "%f", &fin);
			bnew[ibnew].y2=fin;
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].x3=fin;
            fscanf_s(fp, "%f", &fin);
			bnew[ibnew].y3=fin;
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].x4=fin;
            fscanf_s(fp, "%f", &fin);
			bnew[ibnew].y4=fin;

			// свойства материала.
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].rho=fin;
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].cp=fin;
			fscanf_s(fp, "%f", &fin);
			bnew[ibnew].lambda=fin;

			ibnew++;


		}

	}
		

	fclose(fp); // закрытие файла
	return 0;
}

// выясняет в каком блоке находится точка (xpos, ypos).
void inregblock(int lb, BLOCK* b, int lbt, BLOCKU* bt, int lbsec, SECTORU* bsec,
				double xpos, double ypos, double theta, int &it) {
	it=-1;  // изменяется от -1 (не принадлежит блокам) до lb+lbt-1.
	// если лежит от 0 до lb-1 то принадлежит прямоугольному блоку.
	// если лежит от  lb до lbt-1 то принадлежит top блоку.
	for (int iprior=0; iprior<10000; iprior++) {
		// iprior - текущее значение приоритета.
	    for (int i=0; i<lb; i++) {
		    if (b[i].priority==iprior) {
		       if ((xpos>b[i].xS)&&(xpos<b[i].xE)&&(ypos>b[i].yS)&&(ypos<b[i].yE)) {
			      it=i;
		       }
		    }
	    } // прямоугольные блоки.
		for (int i=0; i<lbt; i++) {
			if  (bt[i].priority==iprior) {

				if (bt[i].itoporbottom==CURVETOP) {
				     if (bt[i].icurvetop==MYPARABOLA) {
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>bt[i].yS)&&(ypos<(bt[i].at*xpos*xpos+bt[i].bt*xpos+bt[i].ct))) {
					       it=lb+i;
				        }
				     }
				     else if (bt[i].icurvetop==MYSIN) {
					    const double MPI=3.14159265;
					    if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>bt[i].yS)&&(ypos<(bt[i].dyt+bt[i].epsat*sin(2.0*MPI*bt[i].omegat*(xpos-bt[i].dxt))))) {
					       it=lb+i;
				        }
				     }
					 else if (bt[i].icurvetop==MYCIRCLE) {
						 const double MPI=3.14159265;
						 // внимание если убрать последнее условие проявится скрытая ошибка которая ещё не исправлена.
						 if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>bt[i].yS)&&(ypos<(bt[i].yCt+bt[i].rt*sin(MPI*theta/180.0)))&&(ypos<bt[i].yE)) {
					       it=lb+i;
				        }
					 }
				     else if (bt[i].icurvetop==NONE) {
					    if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>bt[i].yS)&&(ypos<bt[i].yE)) {
			               it=lb+i;
		                }
				     }
				}
				else if (bt[i].itoporbottom==CURVEBOTTOM) {
					if (bt[i].icurvetop==MYPARABOLA) {
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].ab*xpos*xpos+bt[i].bb*xpos+bt[i].cb))&&(ypos<bt[i].yE)) {
					       it=lb+i;
				        }
				     }
				     else if (bt[i].icurvetop==MYSIN) {
					    const double MPI=3.14159265;
					    if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].dyb+bt[i].epsab*sin(2.0*MPI*bt[i].omegab*(xpos-bt[i].dxb))))&&(ypos<bt[i].yE)) {
					       it=lb+i;
				        }
				     }
				     else if (bt[i].icurvetop==NONE) {
					    if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>bt[i].yS)&&(ypos<bt[i].yE)) {
			               it=lb+i;
		                }
				     }
				}
				else if (bt[i].itoporbottom==CURVEDOUBLE) {
					if ((bt[i].icurvetop==MYPARABOLA)&&(bt[i].icurvebottom==MYPARABOLA)) {
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].ab*xpos*xpos+bt[i].bb*xpos+bt[i].cb))&&(ypos<(bt[i].at*xpos*xpos+bt[i].bt*xpos+bt[i].ct))) {
					       it=lb+i;
				        }
				     }
					 if ((bt[i].icurvetop==MYPARABOLA)&&(bt[i].icurvebottom==MYSIN)) {
						 const double MPI=3.14159265;
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].dyb+bt[i].epsab*sin(2.0*MPI*bt[i].omegab*(xpos-bt[i].dxb))))&&(ypos<(bt[i].at*xpos*xpos+bt[i].bt*xpos+bt[i].ct))) {
					       it=lb+i;
				        }
				     }
					 if ((bt[i].icurvetop==MYSIN)&&(bt[i].icurvebottom==MYPARABOLA)) {
                        const double MPI=3.14159265;
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].ab*xpos*xpos+bt[i].bb*xpos+bt[i].cb))&&(ypos<(bt[i].dyt+bt[i].epsat*sin(2.0*MPI*bt[i].omegat*(xpos-bt[i].dxt))))) {
					       it=lb+i;
				        }
				     }
					 if ((bt[i].icurvetop==MYSIN)&&(bt[i].icurvebottom==MYSIN)) {
                        const double MPI=3.14159265;
				        if ((xpos>bt[i].xS)&&(xpos<bt[i].xE)&&(ypos>(bt[i].dyb+bt[i].epsab*sin(2.0*MPI*bt[i].omegab*(xpos-bt[i].dxb))))&&(ypos<(bt[i].dyt+bt[i].epsat*sin(2.0*MPI*bt[i].omegat*(xpos-bt[i].dxt))))) {
					       it=lb+i;
				        }
				     }

				}

			}
		}
		for (int i=0; i<lbsec; i++) {
		    if (bsec[i].priority==iprior) {
				const double MPI=3.14159265;
				if (bsec[i].ip==0) {
					// сектор ограниченный двумя окружностями.
				     if ((sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))>bsec[i].rmin)&&
					    (sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))<bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax))/*
					    (180.0*atan((ypos-bsec[i].yC)/(xpos-bsec[i].xC))/MPI>bsec[i].thetamin)&&
					    (180.0*atan((ypos-bsec[i].yC)/(xpos-bsec[i].xC))/MPI<bsec[i].thetamax))*/ {
			            it=lb+lbt+i;
		             }
				}
				if (bsec[i].ip==1) {
					// внутренний радиус заменяется горизотальной или вертикальной прямой.
					if ((theta>-45.0)&&(theta<45.0)) {
						// внутренний радиус ограничен вертикальной линией.
                        if ((xpos > bsec[i].rmin)&&
					        (sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))<bsec[i].rmax)&&
					        (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax))
					        {
			                   it=lb+lbt+i;
		                    }
					}
					else if ((theta>-135.0)&&(theta<-45.0)) {
						// внутренний радиус ограничен горизонтальной линией
						if ((ypos<bsec[i].rmin)&&
					        (sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))<bsec[i].rmax)&&
					        (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax))
					        {
			                   it=lb+lbt+i;
		                    }
					}
					else if ((theta>45.0)&&(theta<135.0)) {
						// внутренний радиус ограничен горизонтальной линией
						if ((ypos>bsec[i].rmin)&&
					        (sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))<bsec[i].rmax)&&
					        (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax))
					        {
			                   it=lb+lbt+i;
		                    }
					}
					else if ((theta>135.0)&&(theta<225.0)) {
						// внутренний радиус ограничен вертикальной линией.
						if ((xpos<bsec[i].rmin)&&
					        (sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))<bsec[i].rmax)&&
					        (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax))
					        {
			                   it=lb+lbt+i;
		                    }
					}
				}
				if (bsec[i].ip==2) {
					// внешний радиус заменяется горизотальной или вертикальной прямой.
					if ((theta>-45.0)&&(theta<45.0)) {
						// внешний радиус ограничен вертикальной линией.
						if ((sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))>bsec[i].rmin)&&
					    (xpos<bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
					else if ((theta>-135.0)&&(theta<-45.0)) {
						// внешний радиус ограничен горизонтальной линией
						if ((sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))>bsec[i].rmin)&&
					    (ypos>bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
					else if ((theta>45.0)&&(theta<135.0)) {
						// внешний радиус ограничен горизонтальной линией
						if ((sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))>bsec[i].rmin)&&
					    (ypos<bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
					else if ((theta>135.0)&&(theta<225.0)) {
						// внешний радиус ограничен вертикальной линией.
						if ((sqrt((xpos-bsec[i].xC)*(xpos-bsec[i].xC)+(ypos-bsec[i].yC)*(ypos-bsec[i].yC))>bsec[i].rmin)&&
					    (xpos>bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
				}
				if (bsec[i].ip==3) {
					// внутренний и внешний радиус заменяется горизотальной или вертикальной прямой.
                    if ((theta>-45.0)&&(theta<45.0)) {
                        // внешний и внутренний радиусы ограничены вертикальной линией.
						if ((xpos > bsec[i].rmin)&&
					       (xpos<bsec[i].rmax)&&
					       (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                   it=lb+lbt+i;
		                   }
					}
					else if ((theta>-135.0)&&(theta<-45.0)) {
						// внешний и внутренний  радиусы ограничены горизонтальной линией
						if ((ypos<bsec[i].rmin)&&
					    (ypos>bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
					else if ((theta>45.0)&&(theta<135.0)) {
						// внешний и внутренний радиусы ограничен горизонтальной линией
						if ((ypos>bsec[i].rmin)&&
					    (ypos<bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
					else if ((theta>135.0)&&(theta<225.0)) {
						// внешний и внутренний радиусы ограничен вертикальной линией.
						if ((xpos<bsec[i].rmin)&&
					    (xpos>bsec[i].rmax)&&
					    (theta>bsec[i].thetamin)&&(theta<=bsec[i].thetamax)) {
			                 it=lb+lbt+i;
		                }
					}
				}
		    }
	    } // прямоугольные блоки.
	}
	// it номер блока в котором находится данная точка.
} // inregblock

// добавляет новую вершину ели её ещё не было.
void addvertexnumber(int* &unicvertexid, int* &ix, int* &iy, int isize, int ivertex, int ip, int jp) {
    bool found=false;

	for (int i=0; i<isize; i++) {
		if (unicvertexid[i]==ivertex) found=true;
	}

	if (!found) {
		int i=1;
		while (unicvertexid[i]>-1) {
			i++;
		}
		unicvertexid[i]=ivertex;
		ix[i]=ip;
		iy[i]=jp;
	}
}

// добавляет новую вершину если её ещё не было.
// структурированная скошенная четырёхугольная сетка.
void addvertexnumberreal(int* &unicvertexid, Real* &rx, Real* &ry, int isize, int ivertex, Real rpx, Real rpy, Real myepsilon) {
    bool found=false;
    //const double myepsilon=0.0001; //1e-20; // точность определения вещественного нуля.
	// ivertex - любое положительное число, оно по сути нужно только для того чтобы распознать конец массива.


	for (int i=0; i<isize; i++) {
		/*if (unicvertexid[i]==ivertex) found=true;*/
		if ((fabs(rx[i]-rpx)<myepsilon)&&(fabs(ry[i]-rpy)<myepsilon)) found=true;
	}

	if (!found) {
		int i=1;
		while (unicvertexid[i]>-1) {
			i++;
		}
		unicvertexid[i]=ivertex;
		rx[i]=rpx;
		ry[i]=rpy;
	}
} // addvertexnumberreal

// ищет порядковый номер добавленой вершины в списке уникальных вершин.
int findvertex(int* &unicvertexid,int isize, int ivertex) {
	int ir=-1; // если возвратит это значение то такойвершины просто несуществует.
	for (int i=0; i<isize; i++) {
		if (unicvertexid[i]==ivertex) {
			ir=i;
			break; // досрочное прерывание цикла.
		}
	}
	return ir;
	// нумерация всегда начинается с 1 т.к. нулевой индекс всегда заполнен -1.
}

// ищет порядковый номер добавленой вершины в списке уникальных вершин.
int findvertexreal(int* &unicvertexid, Real* &rx, Real* &ry, int isize, Real rpx, Real rpy, Real myepsilon) {
	int ir=-1; // если возвратит это значение то такойвершины просто несуществует.

	//const double myepsilon=0.0001; //1e-20; // точность определения вещественного нуля.

	for (int i=0; i<isize; i++) {
		/*
		if (unicvertexid[i]==ivertex) {
			ir=i;
			break; // досрочное прерывание цикла.
		}*/
        if ((fabs(rx[i]-rpx)<myepsilon)&&(fabs(ry[i]-rpy)<myepsilon)) {
			ir=i;
			break; // досрочное прерывание цикла.
		}
	}
	return ir;
	// нумерация всегда начинается с 1 т.к. нулевой индекс всегда заполнен -1.
}

// постронение структуры сетки для прямоугольных элементов.
void constructtaskdat(MYTASK_DATA &taskdat,
					  Real* xpos, Real* ypos,
					  int inx, int iny, 
					  int lb, BLOCK* b, 
					  int lbt, BLOCKU* bt,
					  int lbsec, SECTORU* bsec) {


	taskdat.nve=4; // четырёхугольные элементы.
	// отображение (i,j)->k  : k=i+j*(inx+1).
	// i=0..inx, j=0..iny.

	int isizeuvid=(inx+1)*(iny+1)+1;
	int* unicvertexid=new int [isizeuvid];
	int* ix=new int [isizeuvid];
	int* iy=new int [isizeuvid];
	
	int **nvtx=new int*[taskdat.nve];
	for (int i=0; i<taskdat.nve; i++) {
		nvtx[i]=new int[isizeuvid];
	}

	Real *rho = new Real [isizeuvid];
	Real *cp =  new Real [isizeuvid];
	Real *lam = new Real [isizeuvid];

	for (int i=0; i<isizeuvid; i++) {
		unicvertexid[i]=-1; // несуществующие номера вершин.

		// несуществующие свойства материалов.
		rho[i]=-1.0;
		cp[i]=-1.0;
		lam[i]=-1.0;
	}

	int ie=1;

	// сканируем сетку :
	for (int i=0; i<inx; i++) {
		for (int j=0; j<iny; j++) {
			Real xc, yc; // точка пересечения диагоналей.
			xc=0.5*(xpos[i]+xpos[i+1]);
			yc=0.5*(ypos[j]+ypos[j+1]);
			int it=-1;
			inregblock(lb, b, lbt, bt, lbsec, bsec, xc, yc, -10000.0, it);
			if ((it>=0)&&(it<lb)) {
				if (b[it].inhollow==1) {
					// не hollow блок.

			        // против часовой стрелки упорядочены.
			        addvertexnumber(unicvertexid, ix, iy, isizeuvid, i+j*(inx+1), i, j);
			        addvertexnumber(unicvertexid, ix, iy, isizeuvid, i+1+j*(inx+1), i+1, j);
			        addvertexnumber(unicvertexid, ix, iy, isizeuvid, i+1+(j+1)*(inx+1), i+1, j+1);
			        addvertexnumber(unicvertexid, ix, iy, isizeuvid, i+(j+1)*(inx+1), i, j+1);
 
				    nvtx[0][ie]=findvertex(unicvertexid,isizeuvid, i+j*(inx+1)); // индексация начинается с единицы (это важно)
				    nvtx[1][ie]=findvertex(unicvertexid,isizeuvid, i+1+j*(inx+1));
				    nvtx[2][ie]=findvertex(unicvertexid,isizeuvid, i+1+(j+1)*(inx+1));
				    nvtx[3][ie]=findvertex(unicvertexid,isizeuvid, i+(j+1)*(inx+1));

				    rho[ie]=b[it].rho;
				    cp[ie]=b[it].cp;
				    lam[ie]=b[it].lambda;

				    ie++;
				}
			}


		}
	}

	// заполнение координат узлов:
	// начало.
	int ic1=0;
	for (int i=0; i<isizeuvid; i++) {
		if (unicvertexid[i]>-1) {
			ic1++;
		}
	}
	taskdat.nodes=ic1;
	taskdat.maxnod=ic1;
	taskdat.x=new Real[ic1];
	taskdat.y=new Real[ic1];
	for (int i=0; i<ic1; i++) {
		taskdat.x[i]=xpos[ix[i+1]];
		taskdat.y[i]=ypos[iy[i+1]];
	}
	// конец.

	// заполнение элементов.
	// начало.
	taskdat.nelmts=ie-1;
	taskdat.maxelm=ie-1;
	taskdat.rho=new Real[ie-1];
	taskdat.cp=new Real[ie-1];
	taskdat.lam=new Real[ie-1];
	taskdat.nvtx=new int*[taskdat.nve];
	for (int i=0; i<taskdat.nve; i++) {
		taskdat.nvtx[i]=new int[ie-1];
	}

	for (int i=0; i<ie-1; i++) {
		taskdat.rho[i]=rho[i+1];
		taskdat.cp[i]=cp[i+1];
		taskdat.lam[i]=lam[i+1];

		taskdat.nvtx[0][i]=nvtx[0][i+1];
		taskdat.nvtx[1][i]=nvtx[1][i+1];
		taskdat.nvtx[2][i]=nvtx[2][i+1];
		taskdat.nvtx[3][i]=nvtx[3][i+1];
	}
	// конец.

	// Грниные условия в области по умолчанию.
	// начало.

	// Граничные условия по умолчанию.
	taskdat.constr=new bool*[3];
	taskdat.potent=new Real*[3];
	for (int i=0; i<VAR_COUNT; i++) {
		taskdat.constr[i]=new bool[ic1];
		taskdat.potent[i]=new Real[ic1];
	}

	// инициализация. Всюду однородные условия Неймана.
	for (int i=0; i<VAR_COUNT; i++) {
		for (int j=0; j<ic1; j++) {
			taskdat.constr[i][j]=false; // всюду условия Неймана.
			taskdat.potent[i][j]=0.0;
		}
	}

	// конец.


	// Освобождение оперативной памяти.
	delete ix;
	delete iy;
	delete unicvertexid;

	delete rho;
	delete cp;
	delete lam;

	for (int i=0; i<taskdat.nve; i++) {
		delete nvtx[i];
	}
	delete nvtx;

} // constructtaskdat

// постронение структуры сетки для top элемента.
void constructtaskdattopelem(MYTASK_DATA &taskdat, 
					  int lb, BLOCK* b, 
					  int lbt, BLOCKU* bt,
					  int lbsec, SECTORU* bsec,
					  int lbnew, BLOCKNEW* bnew,
					   Real* xpos) 
{
	taskdat.nve=4; // четырёхугольные элементы.
	// отображение (i,j)->k  : k=i+j*(inx+1).
	// i=0..inx, j=0..iny.

	int inxgl=300;
	int inygl=40;

	int isizeuvid=(inxgl+1)*(inygl+1)+1;
    int* unicvertexid=new int [isizeuvid];
	double* rx=new double [isizeuvid];
	double* ry=new double [isizeuvid];
	for (int i=0; i<isizeuvid; i++) {
		rx[i]=1e20; // несуществующие значения
		ry[i]=1e20; // инициализация.
	}
	
	int **nvtx=new int*[taskdat.nve];
	for (int i=0; i<taskdat.nve; i++) {
		nvtx[i]=new int[isizeuvid];
	}

	Real *rho = new Real [isizeuvid];
	Real *cp =  new Real [isizeuvid];
	Real *lam = new Real [isizeuvid];

	for (int i=0; i<isizeuvid; i++) {
		unicvertexid[i]=-1; // несуществующие номера вершин.

		// несуществующие свойства материалов.
		rho[i]=-1.0;
		cp[i]=-1.0;
		lam[i]=-1.0;
	}

	int ie=1;

	for (int idnew=0; idnew<lbnew; idnew++) {
		int inx=bnew[idnew].n1;
		int iny=bnew[idnew].n2;
		

		// сканируем сетку :
	    for (int i=0; i<inx; i++) {
		    for (int j=0; j<iny; j++) {

				// Сначала разбиваем эталонный квадрат (-1 +1)x(-1 +1)

				Real x1iso, y1iso, x2iso, y2iso, x3iso, y3iso, x4iso, y4iso;
                Real x1, y1, x2, y2, x3, y3, x4, y4;
			    
				x1iso=-1.0+i*((double)(2.0/inx));
				x2iso=-1.0+(i+1)*((double)(2.0/inx));
				x3iso=-1.0+(i+1)*((double)(2.0/inx));
				x4iso=-1.0+i*((double)(2.0/inx));

				y1iso=-1.0+j*((double)(2.0/iny));
				y2iso=-1.0+(j)*((double)(2.0/iny));
				y3iso=-1.0+(j+1)*((double)(2.0/iny));
				y4iso=-1.0+(j+1)*((double)(2.0/iny));

				Real xA, xB, xC, xD, yA, yB, yC, yD;
				xA=bnew[idnew].x1;
				xB=bnew[idnew].x2;
				xC=bnew[idnew].x3;
				xD=bnew[idnew].x4;

				yA=bnew[idnew].y1;
				yB=bnew[idnew].y2;
				yC=bnew[idnew].y3;
				yD=bnew[idnew].y4;

				x1=xA*0.25*(1.0+x1iso)*(1.0+y1iso)+xB*0.25*(1.0-x1iso)*(1.0+y1iso)+xC*0.25*(1.0-x1iso)*(1.0-y1iso)+xD*0.25*(1.0+x1iso)*(1.0-y1iso);
				x2=xA*0.25*(1.0+x2iso)*(1.0+y2iso)+xB*0.25*(1.0-x2iso)*(1.0+y2iso)+xC*0.25*(1.0-x2iso)*(1.0-y2iso)+xD*0.25*(1.0+x2iso)*(1.0-y2iso);
				x3=xA*0.25*(1.0+x3iso)*(1.0+y3iso)+xB*0.25*(1.0-x3iso)*(1.0+y3iso)+xC*0.25*(1.0-x3iso)*(1.0-y3iso)+xD*0.25*(1.0+x3iso)*(1.0-y3iso);
				x4=xA*0.25*(1.0+x4iso)*(1.0+y4iso)+xB*0.25*(1.0-x4iso)*(1.0+y4iso)+xC*0.25*(1.0-x4iso)*(1.0-y4iso)+xD*0.25*(1.0+x4iso)*(1.0-y4iso);
			    
				y1=yA*0.25*(1.0+x1iso)*(1.0+y1iso)+yB*0.25*(1.0-x1iso)*(1.0+y1iso)+yC*0.25*(1.0-x1iso)*(1.0-y1iso)+yD*0.25*(1.0+x1iso)*(1.0-y1iso);
				y2=yA*0.25*(1.0+x2iso)*(1.0+y2iso)+yB*0.25*(1.0-x2iso)*(1.0+y2iso)+yC*0.25*(1.0-x2iso)*(1.0-y2iso)+yD*0.25*(1.0+x2iso)*(1.0-y2iso);
				y3=yA*0.25*(1.0+x3iso)*(1.0+y3iso)+yB*0.25*(1.0-x3iso)*(1.0+y3iso)+yC*0.25*(1.0-x3iso)*(1.0-y3iso)+yD*0.25*(1.0+x3iso)*(1.0-y3iso);
				y4=yA*0.25*(1.0+x4iso)*(1.0+y4iso)+yB*0.25*(1.0-x4iso)*(1.0+y4iso)+yC*0.25*(1.0-x4iso)*(1.0-y4iso)+yD*0.25*(1.0+x4iso)*(1.0-y4iso);
			    
                 
				Real xc, yc; // центр масс выпуклого четырёхугольника.
                xc=0.25*(x1+x2+x3+x4);
	            yc=0.25*(y1+y2+y3+y4);

				if (bnew[idnew].inhollow==1) {
					// не hollow блок.

					double myeps=0.1*(fmin(fmin(fmin(fabs(x1-xc),fabs(x2-xc)),fmin(fabs(x3-xc),fabs(x4-xc))),fmin(fmin(fabs(y1-yc),fabs(y2-yc)),fmin(fabs(y3-yc),fabs(y4-yc)))));

			        // против часовой стрелки упорядочены.
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+j*(inx+1), x1, y1, myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+j*(inx+1), x2, y2, myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+(j+1)*(inx+1), x3, y3,  myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+(j+1)*(inx+1), x4, y4,  myeps);
 
				    nvtx[0][ie]=findvertexreal(unicvertexid, rx, ry, isizeuvid, x1,y1, myeps); // индексация начинается с единицы (это важно)
				    nvtx[1][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x2,y2, myeps);
				    nvtx[2][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x3,y3, myeps);
				    nvtx[3][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x4,y4, myeps);


				    rho[ie]=bnew[idnew].rho;
				    cp[ie]=bnew[idnew].cp;
				    lam[ie]=bnew[idnew].lambda;

				    ie++;
				}

		    }
	    }




	}


	for (int idtop=0; idtop<lbt; idtop++) {

	int inx=bt[idtop].inx;
	int iny=bt[idtop].iny;

	// сканируем сетку :
	for (int i=0; i<inx; i++) {
		for (int j=0; j<iny; j++) {

			Real x1, y1, x2, y2, x3, y3, x4, y4;
			if (xpos==NULL) {
				x1=bt[idtop].xS+i*((double)((bt[idtop].xE-bt[idtop].xS)/inx));
				x2=bt[idtop].xS+(i+1)*((double)((bt[idtop].xE-bt[idtop].xS)/inx));
				x3=bt[idtop].xS+(i+1)*((double)((bt[idtop].xE-bt[idtop].xS)/inx));
				x4=bt[idtop].xS+i*((double)((bt[idtop].xE-bt[idtop].xS)/inx));
			}
			else {
			    x1=xpos[i];
			    x2=xpos[i+1];
			    x3=xpos[i+1];
			    x4=xpos[i];
			}
			if (bt[idtop].itoporbottom==NONE) {
				y1=bt[idtop].yS+j*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				y2=bt[idtop].yS+(j)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				y3=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				y4=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
			}
			else if (bt[idtop].itoporbottom==CURVETOP) {
			    if (bt[idtop].icurvetop==MYPARABOLA) {
			        y1=bt[idtop].yS+j*((double)(((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-bt[idtop].yS)/iny)));
			        y2=bt[idtop].yS+j*((double)(((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-bt[idtop].yS)/iny)));
			        y3=bt[idtop].yS+(j+1)*((double)(((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-bt[idtop].yS)/iny)));
			        y4=bt[idtop].yS+(j+1)*((double)(((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-bt[idtop].yS)/iny)));
			    }
			    else if (bt[idtop].icurvetop==MYSIN) {
				    const double MPI=3.14159265;
				    y1=bt[idtop].yS+j*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-bt[idtop].yS)/iny)));
				    y2=bt[idtop].yS+j*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-bt[idtop].yS)/iny)));
				    y3=bt[idtop].yS+(j+1)*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-bt[idtop].yS)/iny)));
				    y4=bt[idtop].yS+(j+1)*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-bt[idtop].yS)/iny)));
			    }
				else if (bt[idtop].icurvetop==MYCIRCLE) {
					const double MPI=3.14159265;
					y1=bt[idtop].yS+j*((double)(((bt[idtop].yCt+bt[idtop].rt*sin(MPI*(bt[idtop].angleSt+i*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx)))/180.0)-bt[idtop].yS)/iny)));
				    y2=bt[idtop].yS+j*((double)(((bt[idtop].yCt+bt[idtop].rt*sin(MPI*(bt[idtop].angleSt+(i+1)*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx)))/180.0)-bt[idtop].yS)/iny)));
				    y3=bt[idtop].yS+(j+1)*((double)(((bt[idtop].yCt+bt[idtop].rt*sin(MPI*(bt[idtop].angleSt+(i+1)*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx)))/180.0)-bt[idtop].yS)/iny)));
				    y4=bt[idtop].yS+(j+1)*((double)(((bt[idtop].yCt+bt[idtop].rt*sin(MPI*(bt[idtop].angleSt+i*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx)))/180.0)-bt[idtop].yS)/iny)));
				}
			    else if (bt[idtop].icurvetop==NONE) {
				    y1=bt[idtop].yS+j*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y2=bt[idtop].yS+(j)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y3=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y4=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
			    }
			}
			else if (bt[idtop].itoporbottom==CURVEBOTTOM) {
			    if (bt[idtop].icurvebottom==MYPARABOLA) {
			        y1=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+j*((double)(((bt[idtop].yE-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			        y2=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+j*((double)(((bt[idtop].yE-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y3=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].yE-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y4=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].yE-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			    }
			    else if (bt[idtop].icurvebottom==MYSIN) {
				    const double MPI=3.14159265;
				    y1=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+j*((double)((bt[idtop].yE-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
				    y2=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+j*((double)((bt[idtop].yE-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y3=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].yE-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y4=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].yE-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
			    }
			    else if (bt[idtop].icurvebottom==NONE) {
				    y1=bt[idtop].yS+j*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y2=bt[idtop].yS+(j)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y3=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
				    y4=bt[idtop].yS+(j+1)*((double)((bt[idtop].yE-bt[idtop].yS)/iny));
			    }
			}
			else if (bt[idtop].itoporbottom==CURVEDOUBLE) {
                 if ((bt[idtop].icurvebottom==MYPARABOLA)&&(bt[idtop].icurvetop==MYPARABOLA)) {
			        y1=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+j*((double)(((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			        y2=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+j*((double)(((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y3=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y4=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			     }
				 else if ((bt[idtop].icurvebottom==MYSIN)&&(bt[idtop].icurvetop==MYSIN)) {
				    const double MPI=3.14159265;
				    y1=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+j*((double)((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
				    y2=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+j*((double)((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y3=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y4=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
			    }
				 if ((bt[idtop].icurvebottom==MYPARABOLA)&&(bt[idtop].icurvetop==MYSIN)) {
					const double MPI=3.14159265;
			        y1=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+j*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			        y2=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+j*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y3=(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x2-bt[idtop].dxt))-(bt[idtop].ab*x2*x2+bt[idtop].bb*x2+bt[idtop].cb))/iny)));
			        y4=(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb)+(j+1)*((double)(((bt[idtop].dyt+bt[idtop].epsat*sin(2.0*MPI*bt[idtop].omegat*(x1-bt[idtop].dxt))-(bt[idtop].ab*x1*x1+bt[idtop].bb*x1+bt[idtop].cb))/iny)));
			    }
				else if ((bt[idtop].icurvebottom==MYSIN)&&(bt[idtop].icurvetop==MYPARABOLA)) {
				    const double MPI=3.14159265;
				    y1=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+j*((double)((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
				    y2=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+j*((double)((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y3=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].at*x2*x2+bt[idtop].bt*x2+bt[idtop].ct-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x2-bt[idtop].dxb))))/iny));
				    y4=(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb)))+(j+1)*((double)((bt[idtop].at*x1*x1+bt[idtop].bt*x1+bt[idtop].ct-(bt[idtop].dyb+bt[idtop].epsab*sin(2.0*MPI*bt[idtop].omegab*(x1-bt[idtop].dxb))))/iny));
			    }
			}

			Real xc, yc; // центр масс выпуклого четырёхугольника.
            xc=0.25*(x1+x2+x3+x4);
	        yc=0.25*(y1+y2+y3+y4);
			
			int it=-1;
			if (bt[idtop].icurvetop==MYCIRCLE) {
				// здесь важен угол
				inregblock(lb, b, lbt, bt, lbsec, bsec, xc, yc,bt[idtop].angleSt+0.5*(i*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx)) +(i+1)*((double)((bt[idtop].angleEt-bt[idtop].angleSt)/inx))), it);
			}
			else {
			   inregblock(lb, b, lbt, bt, lbsec, bsec, xc, yc, -10000.0, it);
			}
			if (it==idtop+lb) {
				if (bt[idtop].inhollow==1) {
					// не hollow блок.

					double myeps=0.1*(fmin(fmin(fmin(fabs(x1-xc),fabs(x2-xc)),fmin(fabs(x3-xc),fabs(x4-xc))),fmin(fmin(fabs(y1-yc),fabs(y2-yc)),fmin(fabs(y3-yc),fabs(y4-yc)))));

			        // против часовой стрелки упорядочены.
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+j*(inx+1), x1, y1, myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+j*(inx+1), x2, y2, myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+(j+1)*(inx+1), x3, y3,  myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+(j+1)*(inx+1), x4, y4,  myeps);
 
				    nvtx[0][ie]=findvertexreal(unicvertexid, rx, ry, isizeuvid, x1,y1, myeps); // индексация начинается с единицы (это важно)
				    nvtx[1][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x2,y2, myeps);
				    nvtx[2][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x3,y3, myeps);
				    nvtx[3][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x4,y4, myeps);


				    rho[ie]=bt[idtop].rho;
				    cp[ie]=bt[idtop].cp;
				    lam[ie]=bt[idtop].lambda;

				    ie++;
				}
			}


		}
	}

	}

	for (int idsec=0; idsec<lbsec; idsec++) {

		int inx=bsec[idsec].inr;
		int iny=bsec[idsec].intheta;

	// сканируем сетку :
	for (int i=0; i<inx; i++) { // radius
		for (int j=0; j<iny; j++) { // angle

			Real x1, y1, x2, y2, x3, y3, x4, y4;

			const double MPI=3.14159265;

			double t1, t2;
			t1=MPI*bsec[idsec].thetamin/180.0+j*((double)(MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/iny/180.0));
			t2=MPI*bsec[idsec].thetamin/180.0+(j+1)*((double)(MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/iny/180.0));

			double tavg=bsec[idsec].thetamin+0.5*((j)*((double)((bsec[idsec].thetamax-bsec[idsec].thetamin)/iny))+(j+1)*((double)((bsec[idsec].thetamax-bsec[idsec].thetamin)/iny)));

			if (bsec[idsec].ip==0) {
			
			   x1=bsec[idsec].xC+(bsec[idsec].rmin+i*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*cos(t1);
			   x2=bsec[idsec].xC+(bsec[idsec].rmin+(i+1)*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*cos(t1);
			   x3=bsec[idsec].xC+(bsec[idsec].rmin+(i+1)*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*cos(t2);
			   x4=bsec[idsec].xC+(bsec[idsec].rmin+i*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*cos(t2);
			
			   y1=bsec[idsec].yC+(bsec[idsec].rmin+i*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*sin(t1);
			   y2=bsec[idsec].yC+(bsec[idsec].rmin+(i+1)*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*sin(t1);
			   y3=bsec[idsec].yC+(bsec[idsec].rmin+(i+1)*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*sin(t2);
			   y4=bsec[idsec].yC+(bsec[idsec].rmin+i*((double)((bsec[idsec].rmax-bsec[idsec].rmin)/inx)))*sin(t2);
			}
			else if (bsec[idsec].ip==1) {
				// внутренний радиус прямая.
				
				double rmin1=0.0, rmin2=0.0;

				if ((tavg>-45.0)&&(tavg<45.0)) {
					rmin1=(bsec[idsec].rmin-bsec[idsec].xC)/cos(t1);
					rmin2=(bsec[idsec].rmin-bsec[idsec].xC)/cos(t2);

					x1=bsec[idsec].xC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x2=bsec[idsec].xC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x3=bsec[idsec].xC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);
			        x4=bsec[idsec].xC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);

					rmin1=(bsec[idsec].rmin-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmin2=(bsec[idsec].rmin-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
			
					/*
			        y1=bsec[idsec].yC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
			        y4=bsec[idsec].yC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
					*/

					y1=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        y2=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        y3=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx))));
			        y4=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx))));
				}
				else if ((tavg>-135.0)&&(tavg<-45.0)) {
					//rmin1=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t1);
					//rmin2=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t2);

					rmin1=-(bsec[idsec].rmin-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmin2=-(bsec[idsec].rmin-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);

					x1=bsec[idsec].xC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        x2=bsec[idsec].xC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        x3=bsec[idsec].xC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx))));
			        x4=bsec[idsec].xC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx))));

					rmin1=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t1);
					rmin2=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t2);
			
			        y1=bsec[idsec].yC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
			        y4=bsec[idsec].yC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
				}
				else if ((tavg>45.0)&&(tavg<135.0)) {
					rmin1=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t1);
					rmin2=(bsec[idsec].rmin-bsec[idsec].yC)/sin(t2);

					x1=bsec[idsec].xC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x2=bsec[idsec].xC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x3=bsec[idsec].xC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);
			        x4=bsec[idsec].xC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);
			
			        y1=bsec[idsec].yC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
			        y4=bsec[idsec].yC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
				}
				else if ((tavg>135.0)&&(tavg<225.0)) {
					rmin1=(bsec[idsec].rmin-bsec[idsec].xC)/cos(t1);
					rmin2=(bsec[idsec].rmin-bsec[idsec].xC)/cos(t2);

					x1=bsec[idsec].xC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x2=bsec[idsec].xC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*cos(t1);
			        x3=bsec[idsec].xC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);
			        x4=bsec[idsec].xC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*cos(t2);

                    rmin1=-(bsec[idsec].rmin-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmin2=-(bsec[idsec].rmin-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
			
					/*
			        y1=bsec[idsec].yC+(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*sin(t1);
			        y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
			        y4=bsec[idsec].yC+(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx)))*sin(t2);
					*/

                    y1=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        y2=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((bsec[idsec].rmax-rmin1)/inx))));
			        y3=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((bsec[idsec].rmax-rmin2)/inx))));
			        y4=bsec[idsec].yC+2.0*sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+(i)*((double)((bsec[idsec].rmax-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((bsec[idsec].rmax-rmin2)/inx))));

				}

				

			}
			else if (bsec[idsec].ip==2) {
				// внешний радиус прямая.

				double rmin1=bsec[idsec].rmin, rmin2=bsec[idsec].rmin;
				double rmax1=0.0, rmax2=0.0;


				if ((tavg>-45.0)&&(tavg<45.0)) {

					/*t1=MPI*bsec[idsec].thetamin/180.0+atan(((double)(2.0*j))/((double)(iny)));
					t2=MPI*bsec[idsec].thetamin/180.0+atan(((double)((2.0*j+2.0)))/((double)(iny)));
					tavg=0.5*(t1+t2);
					*/

					

					//rmax1=2.0*(bsec[idsec].rmax-bsec[idsec].xC)/sqrt(2.0);
					//rmax2=2.0*(bsec[idsec].rmax-bsec[idsec].xC)/sqrt(2.0);

					
					rmax1=(bsec[idsec].rmax-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmax2=(bsec[idsec].rmax-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);

					y1=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+i*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((rmax1-rmin1)/inx))));
			        y2=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((rmax1-rmin1)/inx))));
			        y3=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((rmax2-rmin2)/inx))));
			        y4=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+i*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((rmax2-rmin2)/inx))));

					rmax1=(bsec[idsec].rmax-bsec[idsec].xC)/cos(t1);
					rmax2=(bsec[idsec].rmax-bsec[idsec].xC)/cos(t2);

					x1=bsec[idsec].xC+(rmin1+i*((double)((rmax1-rmin1)/inx)))*cos(t1);
			        x2=bsec[idsec].xC+(rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*cos(t1);
			        x3=bsec[idsec].xC+(rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*cos(t2);
			        x4=bsec[idsec].xC+(rmin2+i*((double)((rmax2-rmin2)/inx)))*cos(t2);
				}
				else if ((tavg>-135.0)&&(tavg<-45.0)) {
					rmax1=(bsec[idsec].rmax-bsec[idsec].yC)/sin(t1);
					rmax2=(bsec[idsec].rmax-bsec[idsec].yC)/sin(t2);

                     y1=bsec[idsec].yC+(rmin1+i*((double)((rmax1-rmin1)/inx)))*sin(t1);
			         y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*sin(t1);
			         y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*sin(t2);
			         y4=bsec[idsec].yC+(rmin2+i*((double)((rmax2-rmin2)/inx)))*sin(t2);

                    rmax1=-(bsec[idsec].rmax-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmax2=-(bsec[idsec].rmax-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);



					x1=bsec[idsec].xC+2.0*(cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+i*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i)*((double)((rmax1-rmin1)/inx))));
			        x2=bsec[idsec].xC+2.0*(cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((rmax1-rmin1)/inx))));
			        x3=bsec[idsec].xC+2.0*(cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin1+(i+1)*((double)((rmax1-rmin1)/inx))));
			        x4=bsec[idsec].xC+2.0*(cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+i*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin1+(i)*((double)((rmax1-rmin1)/inx))));

				}
				else if ((tavg>45.0)&&(tavg<135.0)) {
					rmax1=(bsec[idsec].rmax-bsec[idsec].yC)/sin(t1);
					rmax2=(bsec[idsec].rmax-bsec[idsec].yC)/sin(t2);

					y1=bsec[idsec].yC+(rmin1+i*((double)((rmax1-rmin1)/inx)))*sin(t1);
			         y2=bsec[idsec].yC+(rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*sin(t1);
			         y3=bsec[idsec].yC+(rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*sin(t2);
			         y4=bsec[idsec].yC+(rmin2+i*((double)((rmax2-rmin2)/inx)))*sin(t2);


					 rmax1=(bsec[idsec].rmax-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmax2=(bsec[idsec].rmax-bsec[idsec].yC)/sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);

					 x1=bsec[idsec].xC+2.0*cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+i*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((rmax1-rmin1)/inx))));
			        x2=bsec[idsec].xC+2.0*cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((rmax1-rmin1)/inx))));
			        x3=bsec[idsec].xC+2.0*cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((rmax2-rmin2)/inx))));
			        x4=bsec[idsec].xC+2.0*cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0)*((rmin2+i*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((rmax2-rmin2)/inx))));
				}
				else if ((tavg>135.0)&&(tavg<225.0)) {
					//rmax1=-2.0*(bsec[idsec].rmax-bsec[idsec].xC)/sqrt(2.0);
					//rmax2=-2.0*(bsec[idsec].rmax-bsec[idsec].xC)/sqrt(2.0);

					rmax1=-(bsec[idsec].rmax-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);
					rmax2=-(bsec[idsec].rmax-bsec[idsec].xC)/cos(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0);

					 y1=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+i*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+i*((double)((rmax1-rmin1)/inx))));
			         y2=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*(j)/iny-0.5*(rmin1+(i+1)*((double)((rmax1-rmin1)/inx))));
			         y3=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+(i+1)*((double)((rmax2-rmin2)/inx))));
			         y4=bsec[idsec].yC+2.0*(sin(0.5*MPI*(bsec[idsec].thetamax-bsec[idsec].thetamin)/180.0))*((rmin2+i*((double)((rmax2-rmin2)/inx)))*(j+1)/iny-0.5*(rmin2+i*((double)((rmax2-rmin2)/inx))));

					rmax1=(bsec[idsec].rmax-bsec[idsec].xC)/cos(t1);
					rmax2=(bsec[idsec].rmax-bsec[idsec].xC)/cos(t2);

                    x1=bsec[idsec].xC+(rmin1+i*((double)((rmax1-rmin1)/inx)))*cos(t1);
			        x2=bsec[idsec].xC+(rmin1+(i+1)*((double)((rmax1-rmin1)/inx)))*cos(t1);
			        x3=bsec[idsec].xC+(rmin2+(i+1)*((double)((rmax2-rmin2)/inx)))*cos(t2);
			        x4=bsec[idsec].xC+(rmin2+i*((double)((rmax2-rmin2)/inx)))*cos(t2);
				}

				
			
			   

			}
			else if (bsec[idsec].ip==3) {

				// Это ничто иное как трапеция такой элемент должен быть предусмотрен раньше.
				// Конфузор, Диффузор.
				// внешний и внутренний радиус прямая.

				if ((tavg>-45.0)&&(tavg<45.0)) {


				}
				else if ((tavg>-135.0)&&(tavg<-45.0)) {
				}
				else if ((tavg>45.0)&&(tavg<135.0)) {
				}
				else if ((tavg>135.0)&&(tavg<225.0)) {
				}
			}
            
			Real xc, yc; // центр масс выпуклого четырёхугольника.
            xc=0.25*(x1+x2+x3+x4);
	        yc=0.25*(y1+y2+y3+y4);
			
			int it=-1;
			inregblock(lb, b, lbt, bt, lbsec, bsec, xc, yc,tavg, it);
			if (it==idsec+lb+lbt) {
				if (bsec[idsec].inhollow==1) {
					// не hollow блок.

					double myeps=0.1*(fmin(fmin(fmin(fabs(x1-xc),fabs(x2-xc)),fmin(fabs(x3-xc),fabs(x4-xc))),fmin(fmin(fabs(y1-yc),fabs(y2-yc)),fmin(fabs(y3-yc),fabs(y4-yc)))));

					// против часовой стрелки упорядочены.
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+j*(inx+1), x1, y1,myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+j*(inx+1), x2, y2,myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+1+(j+1)*(inx+1), x3, y3,myeps);
			        addvertexnumberreal(unicvertexid, rx, ry, isizeuvid, i+(j+1)*(inx+1), x4, y4,myeps);

					 
				    nvtx[0][ie]=findvertexreal(unicvertexid, rx, ry, isizeuvid, x1,y1,myeps); // индексация начинается с единицы (это важно)
				    nvtx[1][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x2,y2,myeps);
				    nvtx[2][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x3,y3,myeps);
				    nvtx[3][ie]=findvertexreal(unicvertexid, rx, ry,isizeuvid, x4,y4,myeps);


				    rho[ie]=bsec[idsec].rho;
				    cp[ie]=bsec[idsec].cp;
				    lam[ie]=bsec[idsec].lambda;

				    ie++;
				}
			}

		}
	}
	}

	// заполнение координат узлов:
	// начало.
	int ic1=0;
	for (int i=0; i<isizeuvid; i++) {
		if (unicvertexid[i]>-1) {
			ic1++;
		}
	}
	taskdat.nodes=ic1;
	taskdat.maxnod=ic1;
	taskdat.x=new Real[ic1];
	taskdat.y=new Real[ic1];
	for (int i=0; i<ic1; i++) {
		taskdat.x[i]=rx[i+1];
		taskdat.y[i]=ry[i+1];
	}
	// конец.

	// заполнение элементов.
	// начало.
	taskdat.nelmts=ie-1;
	taskdat.maxelm=ie-1;
	taskdat.rho=new Real[ie-1];
	taskdat.cp=new Real[ie-1];
	taskdat.lam=new Real[ie-1];
	taskdat.nvtx=new int*[taskdat.nve];
	for (int i=0; i<taskdat.nve; i++) {
		taskdat.nvtx[i]=new int[ie-1];
	}

	for (int i=0; i<ie-1; i++) {
		taskdat.rho[i]=rho[i+1];
		taskdat.cp[i]=cp[i+1];
		taskdat.lam[i]=lam[i+1];

		taskdat.nvtx[0][i]=nvtx[0][i+1];
		taskdat.nvtx[1][i]=nvtx[1][i+1];
		taskdat.nvtx[2][i]=nvtx[2][i+1];
		taskdat.nvtx[3][i]=nvtx[3][i+1];
	}
	// конец.

	// Грниные условия в области по умолчанию.
	// начало.

	// Граничные условия по умолчанию.
	taskdat.constr=new bool*[3];
	taskdat.potent=new Real*[3];
	for (int i=0; i<VAR_COUNT; i++) {
		taskdat.constr[i]=new bool[ic1];
		taskdat.potent[i]=new Real[ic1];
	}

	// инициализация. Всюду однородные условия Неймана.
	for (int i=0; i<VAR_COUNT; i++) {
		for (int j=0; j<ic1; j++) {
			taskdat.constr[i][j]=false; // всюду условия Неймана.
			taskdat.potent[i][j]=0.0;
		}
	}

	// конец.


	// Освобождение оперативной памяти.
	delete rx;
	delete ry;
	delete unicvertexid;

	delete rho;
	delete cp;
	delete lam;

	for (int i=0; i<taskdat.nve; i++) {
		delete nvtx[i];
	}
	delete nvtx;


} // constructtaskdattopelem

int _tmain(int argc, _TCHAR* argv[])
{

	MYTASK_DATA taskdat;

	int ifine;
	ifine=replace();
	if (ifine!=0) {
		printf("error %d in file obj.txt\n",ifine);
		printf("please, press any key to continue...\n");
		getchar();
	}
	else {

		int inx=33;
		int iny=33;


		double *xpos, *ypos;
		xpos=NULL;
		ypos=NULL;

		int lb, lbt, lbsec, lbnew; // число блоков, универсальных блоков, число секторов.
	    // файл был успешно преобразован.
		readobj(inx, iny, lb, lbt, lbsec, lbnew); // считывание объектов.

		if (lb>0) {
		   // генерирует позиции узлов для серии прямоугольных блоков.
		   simplemeshgen(xpos, ypos, inx, iny, lb, b); 
	       // конструирование прямоугольной сетки.
		   constructtaskdat(taskdat, xpos, ypos, inx, iny, lb, b, lbt, bt, lbsec, bsec);
		}

		if ((lbt>0)||(lbsec>0)||(lbnew>0)) {
		   constructtaskdattopelem(taskdat, lb, b, lbt, bt, lbsec, bsec, lbnew, bnew, NULL);
		}



        exporttecplotxy360(taskdat.nve, taskdat.nelmts, taskdat.nodes, taskdat.nvtx, taskdat.x, taskdat.y, taskdat.potent[TEMP]);


	}

	return 0;
}

