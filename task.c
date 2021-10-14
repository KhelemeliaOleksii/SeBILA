#include<stdio.h>
#include<conio.h>
#include<math.h>
#include <stdlib.h> //malloc
#include "Task.h"
#include <mpi.h>
#include <string.h>	//memcpy

// DEFINITION: Function add inparameter value to outparameter
// IN:	float subsumm - value of partial summ
// OUT: float *sum - value of summ 
void summator(float subsumm, float *sum) {
	*sum = *sum + subsumm;
}
float totalsummator_mpi(int rank, int size, float subSumm) {
	int i;
	float result = 0;
	MPI_Status status;

	if (rank != 0) {
		MPI_Send(&subSumm, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
	} else {
		result += subSumm;
		for(i = 1; i < size; i++) {	
			MPI_Recv(&subSumm, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
			result += subSumm;
		}
	}

	return result;
}

// DEFINITION: Function transforms a set of boundaries
	// to the task format
// IN:	float pointmin_set[] - a set low values of boundaries
	//	float pointmax_set[] - a set up values of boundaries
// OUT: int *count_task - counter of tasks
	//	interval list_task[] - a set  of task
void firsttask_Creator(int dimention, float pointmin_set[], float pointmax_set[], int *count_task, interval list_task[]) {
	int i;
	// add task to list
	list_task[*count_task].key = 0;
	for (i = 0; i < dimention; i++) {
		list_task[*count_task].bou[i][0] = pointmin_set[i];
		list_task[*count_task].bou[i][1] = pointmax_set[i];
	}
	// counter of list +1
	*count_task = *count_task+1;
}
// DEFINITION: Function create first no calculated list of tasks
// IN:	float boundaryX_min[] - a set low values of boundaries
	//	float boundaryX_max[] - a set up values of boundaries
// OUT: int *count_task - counter of tasks
	//	interval list_task[] - a set  of task 
//POST: one can find out structure of the task in Interval_struct.h
	// procedure already have been tested - OK. 
int firsttask_CreatorND_mpi(int rank, int size, int *dimention_key, int dim_iter[], float boundaryX_min[], float boundaryX_max[], int *count_task, interval list_task[]) {
	float intervalX[DIM];	// array of interals in all dimentions
	int i, j;				// iterators
	int task_number1D;		// number of interals in one dimention
	int dim_local;			// depth of nested list

	task_number1D = 3*size;
	dim_local = (DIM-*dimention_key);

	if (dim_local == 1) {
		for (i = 0; i < DIM; i++) {
			intervalX[i] = (boundaryX_max[i] - boundaryX_min[i])/((float)task_number1D); 
		}
		for (i = rank; i < task_number1D; i += size) {
			list_task[*count_task].key = 0;
			dim_iter[dim_local-1] = i;
			for(j=0; j < DIM; j++) {
				list_task[*count_task].bou[j][0] = boundaryX_min[j] + dim_iter[j]*intervalX[j];
				list_task[*count_task].bou[j][1] = boundaryX_min[j] + (dim_iter[j]+1)*intervalX[j];
			}
			*count_task = *count_task+1; 
		}
	} else {
		*dimention_key = *dimention_key + 1;
		for (i = 0; i < task_number1D; i++){
			dim_iter[dim_local-1] = i;
			firsttask_CreatorND_mpi(rank, size, dimention_key, dim_iter, boundaryX_min, boundaryX_max, count_task, list_task);
		}
		*dimention_key = *dimention_key-1;
	}
	return 1;
}


// DEFINITION: Function transforms a set function values
	// to the task format
// USE: reccurtion
// IN:	float funcvalue_set[] - a set of function values	
// OUT: interval list_task[] - a set function values for the task
void one_TaskcreatorND (float funcvalue_set[], int dimention_key, int *next_funccounter, int step,  int granstep_set[], int count_task, interval list_task[]) {
	int granstep_tmp = 0;
	int i;
//	int fcount_tmp;
	if (dimention_key == 0) {
		for (i = 1; i < DIM; i++) {
			granstep_tmp += granstep_set[i];
		}
		for(i = 0; i < 2; i++) {
			list_task[count_task].fvalue[*next_funccounter+i] = funcvalue_set[step + granstep_tmp+i];
		}
		*next_funccounter = *next_funccounter+2;
	} else {
		granstep_set[dimention_key] = 0;
		one_TaskcreatorND(funcvalue_set, dimention_key-1, next_funccounter, step,  granstep_set, count_task, list_task);
		granstep_set[dimention_key] = (int)pow((float)NUMPOINTS1D, dimention_key);
		one_TaskcreatorND(funcvalue_set, dimention_key-1, next_funccounter, step,  granstep_set, count_task, list_task);
	}
}
// DEFINITION: Function transforms sets of boundaries and function values
	// to the task format
// USE: reccurtion
// IN:	float coord_set[][NUMPOINTS1D]	- a set of boundaries
	//	float funcvalue_set[] - a set of function values	
// OUT: int *count_task - counter of tasks
	//	interval list_task[] - a set  of task
// WITHIN: Function uses one_TaskcreatorND() to write function values 
void task_CreatorND(float coord_set[][NUMPOINTS1D], float funcvalue_set[], int *dimention_key, int *step, 
			int dim_iter[], int *count_task, interval list_task[]) {
	int dim_local;	// 1 - the lowest level
					// ...
					// N - the highest level
	int dimstep_correction;	// 0 - the lowest level		0	3	6 /	10	13	16	from the sell 1-4-2-5 to the sell 3-6-4-7 -> step +1
							// 1 - next dimention		1	4	7 /	11	14	17	from the sell 4-7-5-9 to the sell 10-13-11-14 -> step +5
							// 5, 25, 125, 625 ..		2	5	9 /	12	15	19  from ... -> step +25
	int numintervals = NUMPOINTS1D-1; 
	int j, i;
//	int gran_step = 0; 
	int next_funccounter = 0;
	int granstep_set[DIM];  // 0 - the lowest level
							// 1, 5, 25...
							// 5^(N-1) - the highest level
	dim_local = (DIM - *dimention_key);
	if (dim_local == 1) {
		for (i = 0; i < numintervals; i++){
			// a set of iterators for all dimentions
			dim_iter[dim_local-1] = i;
			list_task[*count_task].key = 1;
			// write to task boundaries
			for (j = 0; j <= *dimention_key; j++) {
				list_task[*count_task].bou[j][0] = coord_set[j][dim_iter[j]];
				list_task[*count_task].bou[j][1] = coord_set[j][dim_iter[j]+1];
			}
			// write to task node function value			
			next_funccounter = 0;
			one_TaskcreatorND (funcvalue_set, *dimention_key, &next_funccounter, *step,  granstep_set, *count_task, list_task);
			*step = *step + 1;

			// increase counter of task
			*count_task = *count_task+1;
		}
	} else {
		dimstep_correction = (int)pow((float)NUMPOINTS1D, dim_local-2);
		*dimention_key = *dimention_key + 1;
		for (i = 0; i < numintervals; i++){
			dim_iter[dim_local-1] = i;
			task_CreatorND(coord_set, funcvalue_set, dimention_key, step, dim_iter, count_task, list_task);
			*step = *step + dimstep_correction;
		}
		*dimention_key = *dimention_key-1;
	} 
}

void task_sender(int rank, int size, int series, int *count_task, interval list_task[]) {

	MPI_Status status12r;
	MPI_Request req12r;
	int flag12r = 1; // I'm ready to receive a request for additional tasks
	int countsend;		// count tasks to send
	int k, i, j;		// iterators
	int seriestmp12 = 0;

	MPI_Irecv(&seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);
	MPI_Test(&req12r, &flag12r, &status12r);
	if (!flag12r) {
		MPI_Cancel(&req12r);
	} else {
		if (seriestmp12 == series) {
			MPI_Datatype sndrcvdata, interval_mpi;
			MPI_Datatype oldtypes[4];
			MPI_Aint offsets[4];
			int blklens[4];
			float *bufdatasend;
			int position = 0;
			int intbuf[2];
			interval *sendlist;

			countsend = 0;

			if (*count_task%2) {
				countsend = (*count_task-1)/2;
				*count_task = countsend+1;
				sendlist = (interval*)malloc(countsend*sizeof(interval));
				memcpy(sendlist, &list_task[*count_task+1], countsend*sizeof(interval));
			} else { 
				countsend = (*count_task)/2;																	////
				*count_task = countsend;																		////
				sendlist = (interval*)malloc(countsend*sizeof(interval));
				memcpy(sendlist, &list_task[*count_task+1], countsend*sizeof(interval));
			}
			// error of calculation
			// exclusion of parameter`s value 
			if (countsend == 0) {
				int err = 001;
				fprintf(stdout, "ERROR_%d\n", err);fflush(stdout);
				MPI_Abort(MPI_COMM_WORLD, err);
			}

			// BUILD and SEND preparatory message to receied request:
				//	1) series of the work
				//	2) count of additional tasks
			intbuf[0] = series;			
			intbuf[1] = countsend;		
			MPI_Send(intbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);

			// BUILD and SEND main message to receied request
				//Create a send-recv struct type like Interval structure
			// old mpi data types	
			oldtypes[0] = MPI_INT;									//
			oldtypes[1] = MPI_FLOAT;
			oldtypes[2] = MPI_FLOAT;
			// length of old types
			blklens[0] = 1;														//
			//blklens[1] = pow((float)2, DIM);
			//blklens[2] = pow((float)2, DIM);//
			blklens[1] = 4;
			blklens[2] = 4;//
			// start position of datas
			offsets[0] = 0;
			MPI_Type_size(MPI_INT, &offsets[1]);
			MPI_Type_size(MPI_FLOAT, &offsets[2]);
			offsets[2] = offsets[2]*blklens[1];
			// compound new type
			MPI_Type_struct( 3, blklens, offsets, oldtypes, &interval_mpi );		//
			// declare new mpi type									//
			MPI_Type_commit( &interval_mpi);										//

			// Create a send-recv struct 
			/*Data types*/													//
			oldtypes[0] = MPI_INT;									//
			oldtypes[1] = interval_mpi ;
			// start position of datas
			offsets[0] = 0;
			MPI_Type_size(MPI_INT, &offsets[1]);
			// compound new type
			MPI_Type_struct(2, blklens, offsets, oldtypes, &sndrcvdata);		//
			// declare new mpi type									//
			MPI_Type_commit( &sndrcvdata );										//
			//fprintf(stdout, "Test correct sndrcvdata type create (SENDER)\n");fflush(stdout);

			bufdatasend = (float*)malloc(countsend*sizeof(interval)+sizeof(int));
			position = 0;
			MPI_Pack(&intbuf[0], 1, MPI_INT, bufdatasend, countsend*sizeof(interval)+sizeof(int), &position, MPI_COMM_WORLD);
			MPI_Pack(sendlist, countsend, interval_mpi, bufdatasend, countsend*sizeof(interval)+sizeof(int), &position, MPI_COMM_WORLD);

			MPI_Send(bufdatasend, 1, sndrcvdata, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);
			
			free(bufdatasend);
			free(sendlist);
			MPI_Type_free(&sndrcvdata);
			MPI_Type_free(&interval_mpi);
		}
	}
}

int task_asker(int rank, int size, int series, int *count_task, interval list_task[], int *currentsize, int IncreaseProcId[], int INcountcompare[]){
		// tags for communication
	MPI_Status *status13r;
	MPI_Status status12r;
	MPI_Status status12s;
	MPI_Status status11r;
	MPI_Status status11s;
	MPI_Status status10r;
	MPI_Status status10s;
	MPI_Request *req13r;
	MPI_Request req12s = MPI_REQUEST_NULL;
	MPI_Request req12r;
	MPI_Request req11r = MPI_REQUEST_NULL;
	MPI_Request req11s = MPI_REQUEST_NULL;
	MPI_Request req10r = MPI_REQUEST_NULL;
	MPI_Request req10s = MPI_REQUEST_NULL;
	int *flag13r;
	int flag12s = 1;
	int flag12r = 1; // I'm ready to receive a request for additional tasks
	int flag11r = 1;
	int flag11s = 1;
	int flag10r = 0;
	int flag10s = 1;
	int incount;
	int countsend;		// count tasks to send
	int k, i, j;		// iterators
	int *flagexitIn;	//flag of exit increasing list
//	int *IncreaseProcId;//list of coprocessors Increasing
	int *LeaveProcId;	//list of coprocessors which leave communication
	int LeaveN = 0;		
	int LeaveJ = 0;
	int bye = 1;
	int InN;		
	int InNTemp;
	int Inj;
	int Inprev;
	int flagExitProc = 0;// my flag of exit
	int signexitIncom;	
	int *signallIncom;
	int *signallAsk;
//	int *INcountcompare;	// количество задач в возрастающем потоке
	int seriestmp12 = 0;
	MPI_Datatype sndrcvdata, interval_mpi;
	MPI_Datatype oldtypes[4];
	MPI_Aint offsets[4];
	int blklens[4];
	double *bufdatasend;
	int position = 0;
	int intbuf[2];
	interval *sendlist;

	// list of coprocessors
//	IncreaseProcId = (int*)malloc(size*sizeof(int));
	// list of tasks in processors that we know
	// at the start NULL tasks. We know nothing about.
//	INcountcompare = (int*)malloc(size*sizeof(int));
	// a list of mpi requests for no task
	req13r = (MPI_Request*)malloc(*currentsize*sizeof(MPI_Request));
	// ñòàòóñû ïîëó÷åíûõ ñîîáùåíèé oá îòñóòñòâèè äîïîëíèòåëüíûõ çàäàíèé;
	status13r = (MPI_Status*)malloc(size*sizeof(MPI_Status));
	//a list of flags for exit
	flagexitIn = (int*)malloc(*currentsize*sizeof(int));
	//a list of flags for no additinal task
	flag13r = (int*)malloc(*currentsize*sizeof(int));
	// Ñïèñîê ðàáî÷èõ (èëè çàâåðøèâøèõ ðàáîòó) ïðîöåññîâ;
	LeaveProcId = (int*)malloc(*currentsize*sizeof(int));
	// Ñïèñîê ôëàãîâ íà âõîÿùèå çàäàíèÿ;
	signallIncom = (int *)malloc(*currentsize*sizeof(int));
	// Ñïèñîê ôëàãîâ çàïðîñîâ íà äîïîëíèòå çàäàíèÿ;
	signallAsk = (int *)malloc(*currentsize*sizeof(int));

	// create a list of processors for a conversation
	// myid processor is root one
	// example: size 6, myid 4
	// IncreaseProcId[0] = 4, IncreaseProcId[1] = 5..6..0..1..2..3  
	//for(i = rank; i < size; i++) {
	//	IncreaseProcId[i-rank] = i;
	//}
	//for(i = 0; i < rank; i++) {
	//	IncreaseProcId[size-rank+i] = i;
	//}

	//for (i = 0; i < size; i++) {
	//	fprintf(stdout, "Proc %d coproc %d is %d\n", rank, i, IncreaseProcId[i]);fflush(stdout);
	//}
	InN = *currentsize-1;	// количество элементов в восходящем списке
	InNTemp = InN;	// 
	Inj = 1;		// порядковый номер запрашиваемого процесса(восх),
	Inprev = 0;		// предыдущий запрашиваемый процесс
					//	 в начале им является, собственно, корневой процесс
	
	// list of tasks in processors that we know
	// at the start NULL tasks. We know nothing about.
	// assume all processors have no tasks 
	//for (i = 0; i < size; i++) {
	//	INcountcompare[i] = 0;
	//}
	
	// Null counters:
	// 1) flag list of exit signal
	// 2) flag list of no additianal tasks
	// 3) request list of exit signal

	for (i = 1; i < *currentsize; i++) {																//
		flagexitIn[i] = 0;
		flag13r[i] = 0;
		req13r[i] = MPI_REQUEST_NULL;
	}
	for (k = 1; k < *currentsize; k++) {																//
		signallIncom[k] = 1;																		//
	}

	//
	for (i = 1; i < *currentsize; i++) {																//
		signallAsk[i] = 1;																		//
	}																							//

	fprintf(stdout, "->I'm proc %d. I've got no a task. Current size is %d\n", rank, *currentsize);fflush(stdout);
//		getch();

//---	return 1;
/// е. здесь будут вноситься какие-нибудь изменения, 
/// не забыть их зделать в нижней части программы

//////////////////////////////////////////////////////////////////////////	
	//Выполняется пока получит задания
	// или срабатывает спецыальные прерыватели:
				// 1 - не осталось заданий в соседних процессах
				// 2 - все соседние процессы завершили свою работу
	while(!*count_task) {
		int inbuf[2];
		int outbuf[2]; 
		int byebuf;
		int seriestmp10;

		///
		///проверка сообщений о прекращении работы других процессов
		///
		for(i = 1; i <= InNTemp; i++) {
			if(!flagexitIn[i]){ //запрос на проверку наличия сообщения о выходе от других процессов (один раз)
				MPI_Irecv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				flagexitIn[i] = 1;
			}
			if(!flag13r[i]) { //проверка выполнения запроса на проверку наличия сообщения о выходе от других процессов 
				MPI_Test(&req13r[i], &flag13r[i], &status13r[i]);
			}
			if (flag13r[i]) { // проверка выполнения запроса на проверку наличия сообщения о выходе от других процессов 
								// дала положительный результат
				int seriestmp13;
				seriestmp13 = byebuf;
				fprintf(stdout, "Hello\n");fflush(stdout);

				if (series == seriestmp13) {
//					flag12s = 1;
					//fprintf(stdout, "I'm proc %d. Request 13. Proc %d leave the job\n", rank, status13r[i].MPI_SOURCE);fflush(stdout);
					InNTemp --; //количество соседних процесов уменьшается на 1
					*currentsize = *currentsize-1;
					// е. запрашиваемый на дополнительные задания процес вышел,
					// то запрос на дополнительные задания отменяется
					if (Inj == i) {
						flag11r = 1;
						if (flag10r) {
							MPI_Cancel(&req10r);
							//flag10r = 1;
							flag12s = 1;
						}
					}
					///е. сообщение о прекращении работы процессом принято, тогда
					/// изымаем его из списка рабочих процессов
					if(Inj > i) {
						Inj --;
					}
					for(j = i; j <= InNTemp; j++){
						IncreaseProcId[j] = IncreaseProcId[j+1];
						signallIncom[j] = signallIncom[j+1]; 
					}
					for (j = i; j<=InNTemp; j++) {
						req13r[j] = req13r[j+1];
						flag13r[j] = flag13r[j+1] ;
					}
				} else {
					flagexitIn[i] = 0;
					flag13r[i]=0;
//					req13r[i] = 0;
				}
			}
        }

		///
        /// Сделать вставку на удаление отосланных сообщений процессу, завершившему свою работу
        ///

		/// е. нет работающих процессов, то завершаем работу
        InN = InNTemp; // новое число рабочих процессов;
		if (InN == 0) {//если нет рабочих процессоров, то выход
			//fprintf(stdout, "proc %d leaving caused by no more coprocessors and task\n", rank);	fflush(stdout);
			if (!flag12r) {
				MPI_Cancel(&req12r);
			}
			if (flag10r) {
		        MPI_Cancel(&req10r);
			}
//			fprintf(stdout, "proc %d leaving caused by no more coprocessors and tasks no exist\n", rank);	fflush(stdout);
            return 1;
		}
/////////////////////////////////////////////////

		///
		/// СДЕЛАТЬ и ВЫПОЛНИТЬ запрос на дополнительные задания
		///
		if (flag12s) {
			int tmp = series;
			MPI_Send(&tmp, 1, MPI_INT, IncreaseProcId[Inj], 12, MPI_COMM_WORLD);
//			fprintf(stdout, "I'm proc %d. Send request for add task to proc %d \n", rank, IncreaseProcId[Inj]);fflush(stdout);
			flag12s = 0;
			flag10r = 1;
		}
//		break;
		// получение запрос на дополнительные задания
		// и отсылка сообщения об отсутствии дополнительных заданий
///
/// СДЕЛАТЬ и по-возможности ПОЛУЧИТЬ запрос на дополнительные задания
///
        if (flag12r) {
			MPI_Irecv(&seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);////
            flag12r = 0;																									////
        }
        MPI_Test(&req12r, &flag12r, &status12r);

		// е. запрос на дополнительные задания ПОЛУЧЕН, 
			// тогда СОЗДАТЬ и ПОСЛАТЬ сообщение, что дополнительных заданий НЕТ
			 // ОБНУЛИТЬ счетчик данного процесса
		if (flag12r) {
			if(series == seriestmp12) { // е. серийник полученного сообщения соответствует 
												// серийнику выполняемого общего задания, тогда...
				outbuf[0] = series;
				outbuf[1] = 0;
				MPI_Send(outbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
				countsend = 0;
				flag12r = 1;
				flag10s = 0;
//				fprintf(stdout, "I'm proc %d. No tasks msg send to coproc %d\n", rank, status12r.MPI_SOURCE);fflush(stdout);
			} 
		}
		// СДЕЛАТЬ и по-возможности ПОЛУЧИТЬ сообщение 
			// о КОЛИЧЕСТВЕ запрашиваемых дополнительных заданий
		     //.. получение дополнительных заданий
		if (flag10r) {
			MPI_Irecv(inbuf, 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
			flag10r = 0;
		}
        MPI_Test(&req10r, &flag10r, &status10r);

		/////е. ПОЛУЧЕНО сообщение с количеством дополнительных заданий, 
		///// тогда получить это количество заданий
		/////
        if (flag10r) {
			flag10r = 0;
			flag11r = 1;    // флаг принятия сообщения
					//			MPI_Unpack(intbuf, 2*sizeof(int), &position, &seriestmp, 1, MPI_INT, MPI_COMM_WORLD);
			seriestmp10 = inbuf[0];
			if (series == seriestmp10) 
			{
				incount = inbuf[1];
				/// е. количество равно 0, то обнулить счетчик отправившего процесса
			   if (incount == 0) { // пришло пустое сообщение
					flag12s = 1;
					signallIncom[Inj] = 0; // обозначение отсутствия заданий процесса,
												//отправившего пустое сообщение
					/// изменить порядковый номер запрашиваемого процесса
					/// е. в предыдущем запросе было меньше заданий, пор.номер возрастает
					/// иначе пор.номер спадает
//					fprintf(stdout, "I'm proc %d. Recv no task from %d\n", rank, IncreaseProcId[Inj]);fflush(stdout);

					if ((Inj <= (InN)) && (Inj > 0)) {
						if (incount >= INcountcompare[Inprev] ) {
							Inprev = Inj;
							if (Inj == (InN)) {
								Inj = 1;
							} else {
								Inj++;
							}
						} else {
							Inprev = Inj;
							if (Inj == 1) {
								Inj = (InN);
							} else {
								Inj--;
							}
						}
					}
					INcountcompare[Inprev] = incount;	
			/// е. количество больше 0, то СДЕЛАТЬ и ПОЛУЧИТЬ сообщение с
			/// дополнительными заданиями
			/// активировать счетчик заданий процесса
			   } else if (incount > 0) {
					MPI_Datatype sndrcvdata;
					int seriestmp11;
					int position = 0;
					float *bufdatarecv;
					//////////////////////////////////////////////////////////////////////
					/*//	  Create a send-recv struct type like Interval structure*/				//
					//////////////////////////////////////////////////////////////////////
				   /*Data types*/													//
					// old mpi data types	
					oldtypes[0] = MPI_INT;									//
					oldtypes[1] = MPI_FLOAT;
					oldtypes[2] = MPI_FLOAT;
					// length of old types
					blklens[0] = 1;														//
					//blklens[1] = pow((float)2, DIM);
					//blklens[2] = pow((float)2, DIM);//
					blklens[1] = 4;
					blklens[2] = 4;//
					// start position of datas
					offsets[0] = 0;
					MPI_Type_size(MPI_INT, &offsets[1]);
					MPI_Type_size(MPI_FLOAT, &offsets[2]);
					offsets[2] = offsets[2]*blklens[1];

					// compound new type
					MPI_Type_struct( 3, blklens, offsets, oldtypes, &interval_mpi );		//
					// declare new mpi type									//
					MPI_Type_commit( &interval_mpi);										//

					//////////////////////////////////////////////////////////////////////
					/*//	  Create a send-recv struct */				//
					/////////////////////////////////////////////////////////////////////
					/*Data types*/													//
					oldtypes[0] = MPI_INT;									//
					oldtypes[1] = interval_mpi ;
					// start position of datas
					offsets[0] = 0;
					MPI_Type_size(MPI_INT, &offsets[1]);
					// compound new type
					MPI_Type_struct(2, blklens, offsets, oldtypes, &sndrcvdata);		//
					// declare new mpi type									//
					MPI_Type_commit( &sndrcvdata );										//
					//////////////////////////////////////////////////////////////////////
					bufdatarecv = (float*)malloc(incount*sizeof(interval)+sizeof(int));
					
//					fprintf(stdout, "proc %d recv from %d list of %d(%d) tasks\n", rank, status10r.MPI_SOURCE, incount, inbuf[1]);
					
					MPI_Recv(bufdatarecv, 1, sndrcvdata , IncreaseProcId[Inj], 11, MPI_COMM_WORLD, &status11r);
					MPI_Unpack(bufdatarecv, incount*sizeof(interval)+sizeof(int), &position, &seriestmp11, 1, MPI_INT, MPI_COMM_WORLD);   
					if (series == seriestmp11 ) {
						MPI_Unpack(bufdatarecv, incount*sizeof(interval)+sizeof(int), &position, list_task, incount, interval_mpi, MPI_COMM_WORLD);   
						signallIncom[Inj] = 1;
						*count_task = incount;
						/// изменить порядковый номер запрашиваемого процесса
						/// е. в предыдущем запросе было меньше заданий, пор.номер возрастает
						/// иначе пор.номер спадает
						if ((Inj <= (InN)) && (Inj > 0)) {
							if (incount >= INcountcompare[Inprev] ) {
								Inprev = Inj;
								if (Inj == (InN)) {
									Inj = 1;
								} else {
									Inj++;
								}
							} else {
								Inprev = Inj;
								if (Inj == 1) {
									Inj = (InN);
								} else {
									Inj--;
								}
							}
						}
						INcountcompare[Inprev] = incount;	
					}
					MPI_Type_free(&sndrcvdata);
					free(bufdatarecv);
					if (!flag12r) {
						MPI_Cancel(&req12r);
					}
					fprintf(stdout, "I, %d, have received additional task %d\n", rank, incount);
					return 0;

				} else if (incount < 0) {
					fprintf(stdout, "\nproc %d !!!!!SystemError!!!!\n", rank);
					fflush(stdout);
					MPI_Finalize();
				}
			} 
        }
		// проверить всех ли опросил
        // е. да, то отослать сообщения о своем выходе и выйти
        // иначе продолжить опрашивать
/////
///// е. все счетчики процессов обнулены, то это является достаточным условием 
///// для выхода из программы
///// СДЕЛАТЬ и ВЫПОЛНИТЬ рассылку сообщения о своем выходе
///// ОТМЕНИТЬ запросы на дополнительные задания 
///// и проверку на выход  других процессов
/////
//		//fprintf(stdout, "I'm proc %d. Test flagexit (no cooprocessors have got tasks) %d \n", rank, flagExitProc);fflush(stdout);
//		//getch();
//
        if (!flagExitProc) {
            signexitIncom = 0;
			for (k = 1; k <= InN; k++) {		//
                signexitIncom += signallIncom[k];			//
            }
            if (!signexitIncom) {
                for (i = 1; i <= InNTemp; i++) {
					byebuf = series; 
					MPI_Send(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD);
					MPI_Cancel(&req13r[i]);
				}
				fprintf(stdout, "proc %d leaving caused by no more tasks\n", rank);fflush(stdout);

                // сделать вставку на отмену сделланных запросов(на новые задания, на проверку действующих процессов)
				if (!flag12r) {
                    MPI_Cancel(&req12r);
				}
				//if (!flag10r) {
    //                MPI_Cancel(&req10r);
				//}
                return 1;
//				KEY_exit = 0; // сигнал выхода
//               break;
            }
		}
//////////////////////////////////////////
	}


	free(status13r);
	free(flagexitIn);
	free(req13r);
	free(flag13r);
	free(INcountcompare);
	free(IncreaseProcId);
//	free(list);
	free(signallIncom);
	free(signallAsk);

	fprintf(stdout, "<-Proc %d\n", rank);
	fflush(stdout);


	return 0;
}