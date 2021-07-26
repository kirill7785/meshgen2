// � ����� AliceDataDubl.cpp ���������� ������ � �������
// ��� ������������� � ��������� ������� ��������� AliceFlow.
// ��������� �������� ��������� ��� ���������������� ������ �� ����� ���������� �������������� ������ � ������.
// �������� ��� ����������� ��������� ���������� � ��������� AliceFlow ������ ���� �������� ������.



#define Real double
#define Logical bool

// �� ����� AliceData.h
// begin

#define TEMP 0 // �����������
#define VX 1 // �������������� ��������
#define VY 2 // ������������ ��������
#define VAR_COUNT 3 // ���������� ������� ��� ������� ��������� �������� ������.

typedef struct TMYTASK_DATA {
	int maxnod; // ������������ ����� ���� (����������� �������)
    int maxelm; // ����������� ���������� ����� ��������� (��������� ������ ��� �����)
    // �� ��������� �������� ������������. ��� �������� ������������ ��� ���������� ����� � ������.
    int nve; // ����� ������� ���������� �������� (��������� �������� ������ ��� ����������).

    int **nvtx; // ������ ����� ��� ������� ��������
    // ������������ ��������� ������ �� ������������
    Real *x, *y; // ������� ����������
    //int ***nvtxMG; // ������ ����� ��� ������� �������� ��� ���������� ������� �����������.

    Logical **constr; // ������ ��� ������������� ����������� (���. ����������)
    Real **potent; // ������ ������� �����������
    //Real *rthdsd; // ������ ����� ������� ���������

    // �������� ����������
    Real *rho, *cp, *lam;

	//Real **aelm; // ��� ������ �������� A.
    //Real *rloc; // ��������� ������ �����

    int nodes; // ����� ����� � ������
    int nelmts; // ����� ��������� � ������

} MYTASK_DATA;

// end;

// ���� myexporttecplot.cpp
// ������.

// �������� ���������� �����
// ������� ���������� ������� � ��������� tecplot360
// ������������ �������� � ��� ����������� � ��� ���������� �����.
void exporttecplotxy360(int nve, int nelmts, int nodes, int** nvtx, Real* x, Real* y, Real* potent)
{
	FILE *fp;
	errno_t err;
	// �������� ����� ��� ������.
	if ((err = fopen_s( &fp, "tigel_interior_tec.dat"/*"temp_xy.PLT"*/, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"temperature\"\n");

		// ������ ��� ����������
		//fprintf(fp, "VARIABLES = x, y, \"temperature\" , Vx, Vy, Mag \n");
		if (btriangle) {
			fprintf(fp, "VARIABLES = x, y \n"); // ��� ����������� (����� ������ ������ ��� �� �������� ���������).
		}
		else {
		   fprintf(fp, "VARIABLES = x, y, \"temperature\" \n");
		}

		// ������ ���������� � �����
		if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nodes, nelmts);
        if (nve==4) {
			if (btriangle) {
				// ������ �������������� ��� ��� ������������, ������� 
				// ����������� ��������� � ��� ���� ������.
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nodes, 2*nelmts);
			}
			else {
				// �������� �������������� ��������.
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=QUADRILATERAL, F=FEBLOCK\n\n", nodes, nelmts);
			}
		}

		int i=0; // �������� 
		int j=0; // ����� for

		// ������ x
	    for (i=0; i<nodes; i++) {	
				fprintf(fp, "%e ", x[i]);
				if (i%10==0) fprintf(fp, "\n");
		}
			
		fprintf(fp, "\n");
          
		// ������ y
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", y[i]);
            if (i%10==0) fprintf(fp, "\n");
		}
			
        fprintf(fp, "\n");

		if (!btriangle) {
		    // ������ �����������
		    for (i=0;i<nodes; i++) {
			    //fprintf(fp, "%e ", potent[TEMP][i]);
			    fprintf(fp, "%e ", potent[i]);
                if (i%10==0) fprintf(fp, "\n");
		    }

		    fprintf(fp, "\n");

		}
		/*
		// ������ �������������� ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
        // ������ ������������ ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

		// ������ ������ ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        */

		if (btriangle) {
			// ������ �������������� ������� �� ���� �������������.
			for (i=0;i<nelmts; i++) {
				fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[3][i]);
				fprintf(fp, "%d %d %d\n", nvtx[1][i], nvtx[3][i], nvtx[2][i]);
			}
		}
		else {
		    // ������ ���������� � ���������� �����
		    for (i=0;i<nelmts; i++) {
			    if (nve==3) fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i]);
			    if (nve==4) fprintf(fp, "%d %d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i], nvtx[3][i]);
		    }
		}

		fclose(fp); // �������� �����
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
			// ��������� �������
			FILE *fp1;
	        errno_t err1;
            if ((err1 = fopen_s( &fp1, "gru.txt", "r")) != 0) {
		         printf("Open gru.txt File Error\n");
	        }
	        else {

				FILE *fp2;
	            errno_t err2;
	            // �������� ����� ��� ������.
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
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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
						 // ������ x
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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
						 // ������ y
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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

						 // ���� ������ ����� ���������� ��������� ������� � �������.
						 int iprev=-1;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
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
								 // ������������ �����.
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

// �������� ���������� �����
// ������� ���������� ������� � ��������� tecplot360
// ������������ �������� � ��� ����������� � ��� ���������� �����.
void exporttecplotxy360extrude3D(int nve, int nelmts, int nodes, int** nvtx, Real* x, Real* y, Real* potent)
{
	FILE *fp;
	errno_t err;

	float hz=0.2; // ��� � ������� ���������.
	int nz=10; // ���������� ���� � 3D (�� ����� 2��).

	// �������� ����� ��� ������.
	if ((err = fopen_s( &fp, "tigel_interior_tec.dat"/*"temp_xy.PLT"*/, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		// ������ ���������
		fprintf(fp, "TITLE = \"temperature\"\n");

		// ������ ��� ����������
		//fprintf(fp, "VARIABLES = x, y, \"temperature\" , Vx, Vy, Mag \n");
		if (btriangle) {
			fprintf(fp, "VARIABLES = x, y, z \n"); // ��� ����������� (����� ������ ������ ��� �� �������� ���������).
		}
		else {
		   fprintf(fp, "VARIABLES = x, y, z, \"temperature\" \n");
		}

		// ������ ���������� � �����
		if (nve==3) fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=TRIANGLE, F=FEBLOCK\n\n", nz*nodes, nelmts*(nz-1));
        if (nve==4) {
			if (btriangle) {
				// ������ �������������� ��� ��� ������������, ������� 
				// ����������� ��������� � ��� ���� ������.
				// ���� � ��� ��� �������������� ������ BRICK � TETRAHEDRON
                fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=PRIZM, F=FEBLOCK\n\n", nz*nodes, 2*nelmts*(nz-1));
			}
			else {
				// �������� �������������� ��������.
			    fprintf(fp, "ZONE T=\"Rampant\", N=%d, E=%d, ET=BRICK, F=FEBLOCK\n\n", nz*nodes, nelmts*(nz-1));
			}
		}

		int i=0; // �������� 
		int j=0; // ����� for

		// ������ x
		for (int k=0; k<nz; k++) {
	    for (i=0; i<nodes; i++) {	
				fprintf(fp, "%e ", x[i]);
				if (i%10==0) fprintf(fp, "\n");
		}
		}
			
		fprintf(fp, "\n");
          
		// ������ y
		for (int k=0; k<nz; k++) {
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", y[i]);
            if (i%10==0) fprintf(fp, "\n");
		}
		}
			
        fprintf(fp, "\n");

		// ������ z
		for (int k=0; k<nz; k++) {
		for (i=0;i<nodes; i++) {
		 	fprintf(fp, "%e ", 0.0+k*hz);
            if (i%10==0) fprintf(fp, "\n");
		}
		}

		if (!btriangle) {
		    // ������ �����������
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
		// ������ �������������� ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VX][i]);
            if (i%10==0) fprintf(fp, "\n");
		}

		fprintf(fp, "\n");
        // ������ ������������ ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", potent[VY][i]);
            if (i%10==0) fprintf(fp, "\n");
		}
        fprintf(fp, "\n");

		// ������ ������ ��������
		for (i=0;i<nodes; i++) {
			fprintf(fp, "%e ", sqrt(potent[VX][i]*potent[VX][i]+potent[VY][i]*potent[VY][i]));
            if (i%10==0) fprintf(fp, "\n");
		}			

		fprintf(fp, "\n");
        */

		if (btriangle) {
			// ������ �������������� ������� �� ���� �������������.
			for (int k=0; k<nz-1; k++) {
			for (i=0;i<nelmts; i++) {
				fprintf(fp, "%d %d %d %d %d %d\n", nvtx[0][i]+k*nodes, nvtx[1][i]+k*nodes, nvtx[3][i]+k*nodes,nvtx[0][i]+(k+1)*nodes, nvtx[1][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes);
				fprintf(fp, "%d %d %d %d %d %d\n", nvtx[1][i]+k*nodes, nvtx[3][i]+k*nodes, nvtx[2][i]+k*nodes,nvtx[1][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes, nvtx[2][i]+(k+1)*nodes);
			}
			}
		}
		else {
		    // ������ ���������� � ���������� �����
			for (int k=0; k<nz-1; k++) {
		    for (i=0;i<nelmts; i++) {
			    if (nve==3) fprintf(fp, "%d %d %d\n", nvtx[0][i], nvtx[1][i], nvtx[2][i]);
			    if (nve==4) {
					fprintf(fp, "%d %d %d %d %d %d %d %d\n", nvtx[0][i]+k*nodes, nvtx[1][i]+k*nodes, nvtx[2][i]+k*nodes, nvtx[3][i]+k*nodes, nvtx[0][i]+(k+1)*nodes, nvtx[1][i]+(k+1)*nodes, nvtx[2][i]+(k+1)*nodes, nvtx[3][i]+(k+1)*nodes);
				}
		    }
			}
		}

		fclose(fp); // �������� �����
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
			// ��������� �������
			FILE *fp1;
	        errno_t err1;
            if ((err1 = fopen_s( &fp1, "gru.txt", "r")) != 0) {
		         printf("Open gru.txt File Error\n");
	        }
	        else {

				FILE *fp2;
	            errno_t err2;
	            // �������� ����� ��� ������.
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
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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
						 // ������ x
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", x[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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
						 // ������ y
						 for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
								 if (fabs(y[i]-fiso)<1.0e-20) {
									 if ((x[i]>fstart)&&(x[i]<fend)) {
										 fprintf(fp2, "%e ", y[i]);
				                         if (inod%10==0) fprintf(fp2, "\n");
										 inod++;
									 }
								 }
							 }
							 if (iorient==1) {
								 // ������������ �����.
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

						 // ���� ������ ����� ���������� ��������� ������� � �������.
						 int iprev=-1;
	                     for (i=0; i<nodes; i++) {	
                             if (iorient==0) {
								 // �������������� �����.
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
								 // ������������ �����.
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


// �����


