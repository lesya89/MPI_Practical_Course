// copyright : (C) by lesya89
#include <mpi.h> 
#include <iostream>
#include <ctime>
#include <assert.h>
#include <cstdlib>
#include <cmath>

using namespace std;

void MY_MPI_MINLOC(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				if (((int *)recvbuf)[i] < ((int *)sendbuf)[i])
					((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}
	}
}
void MY_MPI_MAXLOC(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				if (((int *)recvbuf)[i] < ((int *)sendbuf)[i])
					((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}
	}
}
void MY_MPI_LXOR(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] != ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_BXOR(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] ^ ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_BOR(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] | ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_LOR(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] || ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_BAND(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] & ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_LAND(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] && ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_PROD(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] * ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_MIN(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				if (((int *)recvbuf)[i] > ((int *)sendbuf)[i])
					((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}
	}
}
void MY_MPI_MAX(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				if (((int *)recvbuf)[i] < ((int *)sendbuf)[i])
					((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}
	}
}
void MY_MPI_SUMM(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	//std::cout<<"I am in 2 if" << std::endl;
	MPI_Status st;
	int ProcNum;
	MPI_Comm_size(comm, &ProcNum);
	for (int j = 0; j < ProcNum; j++) {
		if (j != root)
			MPI_Recv(sendbuf, count, type, j, 0, comm, &st);
		for (int i = 0; i < count; i++)
		{
			if (j == root) {
				((int *)recvbuf)[i] = ((int *)sendbuf)[i];
				//std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
			}
			else {
				((int *)recvbuf)[i] = ((int *)recvbuf)[i] + ((int *)sendbuf)[i];
				//int tmp = ((int *)recvbuf)[i];
				//std::cout << "result " << j << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i] << std::endl;
			}
		}

	}
}
void MY_MPI_SUMM_Tree(void *sendbuf, void *recvbuf, int count) {
	//std::cout<<"I am in 2 if" << std::endl;
	//MPI_Status st;
	int ProcNum;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcNum);
	for (int i = 0; i < count; i++)
	{
		//if (j == root) {
			//((int *)recvbuf)[i] = ((int *)sendbuf)[i];
			//std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
		//}
		//else {
		((int *)recvbuf)[i] = ((int *)recvbuf)[i] + ((int *)sendbuf)[i];
		//int tmp = ((int *)recvbuf)[i];
	//	std::cout << "result " << count << " " << ((int *)recvbuf)[i] << "  last " << ((int *)sendbuf)[i]<<" rank = " << ProcNum << " i = "<< i << std::endl;
	//}
	}

}


int MY_MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	int ProcNum, ProcRank;
	MPI_Comm_size(comm, &ProcNum);
	MPI_Comm_rank(comm, &ProcRank);
	if (ProcRank != root) {
		//std::cout << "I am in not root" << std::endl;
		MPI_Send(sendbuf, count, type, root, 0, comm);
	}
	else {
		if (ProcRank == root)
		{
			//std::cout << "I am in root " << ProcNum << std::endl;
		//	MPI_Status st;
			//int j = 0;
			if (type == MPI_INT) {
				if (op == MPI_SUM) {
					MY_MPI_SUMM(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_MAX) {
					MY_MPI_MAX(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_MIN) {
					MY_MPI_MIN(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_PROD) {
					MY_MPI_PROD(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_LAND) {
					MY_MPI_LAND(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_LOR) {
					MY_MPI_LOR(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_BAND) {
					MY_MPI_BAND(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_BOR) {
					MY_MPI_BOR(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_BXOR) {
					MY_MPI_BXOR(sendbuf, recvbuf, count, type, op, root, comm);
					return 0;
				}
				if (op == MPI_LXOR) {
					MY_MPI_LXOR(sendbuf, recvbuf, count, type, op, root, comm);
				}
				if (op == MPI_MAXLOC)
					return -1;
				//MY_MPI_MAXLOC(sendbuf, recvbuf, count, type, op, root, comm);
				if (op == MPI_MINLOC)
					return -1;


				//MY_MPI_MINLOC(sendbuf, recvbuf, count, type, op, root, comm);
			}
		}
		

	}
 return -10;
}
void recursion(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm, int* massProcRankSend, int size, int h) {
	int curentSize, procNum, rank;
	int* curentMass;
	int* countMassElement = new int[count];
	MPI_Comm_size(comm, &procNum);
	MPI_Comm_rank(comm, &rank);
	MPI_Status st;
	if (size % 2) {
		curentSize = size / 2 + 1;
		curentMass = new int[curentSize];
	}
	else {
		curentSize = size / 2;
		curentMass = new int[curentSize];
	}
	for (int i = 0; i < curentSize; i++) {
		curentMass[i] = -1;
	}
	h++;
	//int t = 0;
	for (int i = 0, t = 0; i < size; i++) {
		if (massProcRankSend[i] == rank) {
			if (!(size % 2)) {
				if (!(i % 2)) {
					MPI_Send(sendbuf, count, type, massProcRankSend[i + 1], 0, comm);
					for (int q = 0; q < count; q++)
					{
						//std::cout << "result " << count << " " << ((int *)recvbuf)[q] << "  last " << ((int *)sendbuf)[q] << " rank = " << rank << " q = " << q << std::endl;
						countMassElement[q] = ((int *)(recvbuf))[q];
					}
				}
				else {
					MPI_Recv(recvbuf, count, type, massProcRankSend[i - 1], 0, comm, &st);
					MY_MPI_SUMM_Tree(sendbuf, recvbuf, count);
					sendbuf = recvbuf;
				}
			}
			else {

				if ((rank != procNum - 1)) {
					if (!(i % 2)) {
						//std::cout << "i = " << i << " rank = " << rank << " rank send " << massProcRankSend[i + 1] << " curentSize = " << curentSize << std::endl;
						MPI_Send(sendbuf, count, type, massProcRankSend[i + 1], 0, comm);
						for (int q = 0; q < count; q++)
						{
							//std::cout << "result " << count << " " << ((int *)recvbuf)[q] << "  last " << ((int *)sendbuf)[q] << " rank = " << rank << " q = " << q << std::endl;
							countMassElement[q] = ((int *)(recvbuf))[q];
						}
					}
					else {
						for (int r = 0; r < count; r++) {
							((int *)recvbuf)[i] = ((int *)sendbuf)[i];
						}
						MPI_Recv(sendbuf, count, type, massProcRankSend[i - 1], 0, comm, &st);
						MY_MPI_SUMM_Tree(sendbuf, recvbuf, count);
						for (int r = 0; r < count; r++) {
							((int *)sendbuf)[i] = ((int *)recvbuf)[i];
						}

					}
				}
				else {
					curentMass[t] = i;
					t++;
				}

			}
		}
		if (!(size % 2)) {
			if (i % 2) {
				curentMass[t] = i;
				t++;
			}
		}
	}

	if (curentSize != 1) {
		recursion(sendbuf, countMassElement, count, type, op, root, comm, curentMass, curentSize, h);
	}
	else {
		if (curentMass[0] != root) {
			if (rank == curentMass[0])
			{
				MPI_Send(sendbuf, count, type, root, 0, comm);
			}
			if (rank == root) {

				MPI_Recv(recvbuf, count, type, curentMass[0], 0, comm, &st);
			}
		}
		else {
			recvbuf = sendbuf;
		}
	}
}

int MY_MPI_Reduce_Tree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
	int ProcNum, ProcRank;
	MPI_Comm_size(comm, &ProcNum);
	MPI_Comm_rank(comm, &ProcRank);
	int * massProcRankSend = new int[ProcNum];
	
	for (int i = 0; i < ProcNum; i++) {
		massProcRankSend[i] = i;
	}
	recursion(sendbuf, recvbuf, count, type, op, root, comm, massProcRankSend, ProcNum, 0);
	if (ProcRank == root) {

	}
	return 0;
}

int main(int argc, char* argv[]) {
	int n = atoi(argv[1]);
	int *mas = new int[n];
	int *mas_r = new int[n];
	int ProcRank;

	double Time_my1 = 0, Time_my2 = 0,Time_tree1=0, Time_tree2=0, Time_MPI1 = 0, Time_MPI2 = 0;
		
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	std::srand(static_cast<int>(time(NULL)));

	for (int i = 0; i < n; i++) {
		mas[i] = 10 + std::rand() % 1000 + ProcRank;
		//std::cout << " rank = " << rank << " mass[" << i << "] "<<mas[i];
	}
	//std::cout << std::endl;

	if (ProcRank == 0)
	{
		Time_my1 = MPI_Wtime();
	}
	MY_MPI_Reduce(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{

		Time_my2 = MPI_Wtime();
		Time_tree1 = MPI_Wtime();
	}
		
	MY_MPI_Reduce_Tree(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		Time_tree2 = MPI_Wtime();
		//std::cout << " mass resualt ";
		for (int i = 0; i < n; i++) {
			//std::cout << mas_r[i] << " ";
		}
		std::cout << std::endl;
		Time_MPI1 = MPI_Wtime();
	}

	MPI_Reduce(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		Time_MPI2 = MPI_Wtime();
		std::cout << "reduce realisation MPI: t =  " << Time_MPI2 - Time_MPI1 << std::endl;
		std::cout << "reduce realisation my tree: t =  " << Time_tree2 - Time_tree1 << std::endl;
		std::cout << "reduce realisation my: t =  " << Time_my2 - Time_my1 << std::endl;
	}
	MPI_Finalize();
	return 0;
}
