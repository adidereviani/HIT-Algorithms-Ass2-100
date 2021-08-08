/*
Assigned by:
Student1_Liron_Levi #312409592
Student2_Adi_Dereviani #305674731
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <assert.h> //malloc/calloc..
#include <time.h> //random
#include <math.h> //log

#define V 1000 //vertex num

int** build_random_graph(double p);
void freeDynamicMatrix(int **mat, int n); // Dynamic 2D array free function
int diameter(int** graph);
int connectivity(int** graph);
int Is_Isolated(int** graph);

void BFS(int ** graph, int v);

//**Q functions**//
void push_queue(int vertex);
int pop_queue();
int isEmpty_queue();

//**global**//
int visited[V], parent[V], distance[V];
int queue[V], front = -1, rear = -1, dist = 0;

void main()
{
	int i = 0, j = 0;
	int **graph;
	double num_graphs = 500;
	int connected = 0;
	int diam = 0;
	int single_ver = 0;
	double prob_1 = 0, prob_2 = 0, prob_3 = 0;
	double Threshold1 = 0, Threshold2 = 0, Threshold3 = 0; // to print the thresholds
	double count_1 = 0, count_2 = 0, count_3 = 0;
	double not_exist_1[10], not_exist_2[10], not_exist_3[10];  //count the graphs
	double exist_1[10], exist_2[10], exist_3[10]; //count the graphs



	//**FILE**//
	FILE *the_results_1 = fopen("connectivity results.csv", "w+");
	FILE *the_results_2 = fopen("diameter results.csv", "w+");
	FILE *the_results_3 = fopen("Is_Isolated results.csv", "w+");
	if ((the_results_1 == NULL) || (the_results_2 == NULL) || (the_results_3 == NULL))
	{
		perror("Unable to open the files");
		exit(1);
	}


	double p1[10] = { 0.0033275, 0.0041550 , 0.0049825 , 0.0053101 , 0.0066376 , 0.0070376 ,0.0076651 , 0.00802926 , 0.00856201 , 0.0090476 }; //1000 0.006907755279
	double p2[10] = { 0.1138 , 0.1144 , 0.1152 , 0.1167, 0.1172 , 0.1179 , 0.1185 , 0.1190 , 0.125 , 0.135 }; //1000 0.1175394
	double p3[10] = { 0.0046275, 0.0051550 , 0.0058825 , 0.0063101 , 0.0068376 , 0.0073376 ,0.0078651 , 0.00802926 , 0.00876201 , 0.0093476 }; //1000 0.006907755279


	srand(time(NULL)); //random

	//log(V) = ln(V) //show the Thresholds
	Threshold1 = log(V) / V;
	printf("the Threshold1 is : %f\n", Threshold1); //show the connectivity Threshold

	Threshold2 = sqrt((2 * log(V)) / V);
	printf("the Threshold2 is : %f\n", Threshold2); //show the diameter Threshold

	Threshold3 = Threshold1;
	printf("the Threshold3 is : %f\n", Threshold3); //show the Is_Isolated Threshold

	//**attribute 1 connectivity**//
	printf("attribute 1:connectivity\n"); //show starting time
	fprintf(the_results_1, "attribute:  connectivity \n");
	fprintf(the_results_1, " P: , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf \n", p1[0], p1[1], p1[2], p1[3], p1[4], p1[5], p1[6], p1[7], p1[8], p1[9]);
	fprintf(the_results_1, " The probability that the feature connectivity exists: ,");
	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < num_graphs; j++)
		{
			graph = build_random_graph(p1[i]); // return random graph
			connected = connectivity(graph); // check if the graph is connected = 1 / unconnected = 0

			if (connected == 1) //connected 
			{
				count_1 += 1; //count the connected graphs
			}

			freeDynamicMatrix(graph, V);
		}
		prob_1 = count_1 / num_graphs; // (num of graph that connected)  /  (num of graphs that checking)
		not_exist_1[i] = (num_graphs - count_1); //count the graphs that unconnected
		exist_1[i] = count_1; //count the graphs that connected
		fprintf(the_results_1, " %lf ,", prob_1);
		count_1 = 0; // restart count
	}
	printf("connectivity - DONE\n"); // show ending time 
	fclose(the_results_1);


	//**attribute 2 diameter**//
	printf("\n attribute 2: diameter\n"); // show starting time
	fprintf(the_results_2, "attribute:  diameter \n");
	fprintf(the_results_2, " P: , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf \n", p2[0], p2[1], p2[2], p2[3], p2[4], p2[5], p2[6], p2[7], p2[8], p2[9]);
	fprintf(the_results_2, "The probability that the feature diameter <= 2 exists: ,");
	for (i = 0; i < 10; i++) // num of P
	{
		count_2 = 0; //restart counter
		for (j = 0; j < num_graphs; j++)
		{
			graph = build_random_graph(p2[i]);// return random graph
			diam = diameter(graph); //get the diam
			if ((diam <= 2) && (diam != -1)) // -1 its unconnected graph with infinity diam
			{
				count_2 += 1; // count the graphs with diam <= 2
			}
			freeDynamicMatrix(graph, V);
		}
		prob_2 = (count_2 / num_graphs); //the prob to get diam <= 2 for each P
		not_exist_2[i] = (num_graphs - count_2); // count the graphs with diam > 2
		exist_2[i] = count_2;// count the graphs with diam <= 2
		fprintf(the_results_2, " %lf ,", prob_2);
	}
	printf("diameter - DONE\n"); // show ending time
	fclose(the_results_2);

	//**attribute 3 Is_Isolated**//
	printf("\n attribute 3: Is_Isolated\n"); // show starting time
	fprintf(the_results_3, "attribute:  Is_Isolated \n");
	fprintf(the_results_3, " P: , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf \n", p3[0], p3[1], p3[2], p3[3], p3[4], p3[5], p3[6], p3[7], p3[8], p3[9]);
	fprintf(the_results_3, " The probability that the feature Isolated vertex exists: ,");
	for (i = 0; i < 10; i++) // num of P
	{
		count_3 = 0;
		for (j = 0; j < num_graphs; j++)
		{
			graph = build_random_graph(p3[i]);// return random graph
			single_ver = Is_Isolated(graph); // if there is a single vertex = 1 , else = 0
			if (single_ver == 1) //single ver exist
			{
				count_3 += 1; // count the graphs with single vertex
			}
			freeDynamicMatrix(graph, V);
		}
		prob_3 = (count_3 / num_graphs); //the prob that a single ver exist for each P
		not_exist_3[i] = (num_graphs - count_3);// count the graphs that single ver is NOT exist
		exist_3[i] = count_3;// count the graphs that single ver is exist
		fprintf(the_results_3, " %lf ,", prob_3);
	}
	printf("Is_Isolated - DONE\n"); // show ending time

	fclose(the_results_3);
}

int** build_random_graph(double p)
{
	int** graph = NULL;
	int i = 0, j = 0;
	int rows = V, cols = V;
	float temp_p;

	graph = (int**)calloc(rows, sizeof(int));
	assert(graph);
	for (i = 0; i < cols; i++) //dyn matrix
	{
		graph[i] = (int*)calloc(cols, sizeof(int)); //new line
		assert(graph[i]); //check if its work
	}

	for (i = 0; i < cols; i++)
	{
		for (j = i + 1; j < rows; j++)
		{
			temp_p = (double)rand() / RAND_MAX; //Lottery number if its smaller than the p then there is a edge between the vertices
			if ((temp_p <= p) && (i != j))
			{ //update the rows and cols of the vertices (edges)
				graph[j][i] = 1;
				graph[i][j] = 1;
			}
		}
	}
	return (graph);
}


int diameter(int** graph)
{
	int   i = 0, j = 0, k = 0, diam = -1, unconnected = 0;

	for (i = 0; i < V; i++) //check all the vertices
	{
		BFS(graph, i + 1); // call BFS for the i ver
		if (i == 0) //first ver - check uf the graph is unconnected then the diam = infinity (-1)
		{
			for (k = 1; k < V; k++)// check the sidtance from the i ver to each ver in the graph
			{
				if (distance[k] == -1) //if there is ver with infinity diam so the graph is unconnected
				{
					unconnected = 1;
					break;
				}
			}
		}
		if ((dist > diam) && (unconnected == 0)) //update diam
		{
			diam = dist;
		}

		if (unconnected == 1)///diam = infinity (-1) so we stop checking
		{
			break;
		}
	}
	return diam;
}

int Is_Isolated(int** graph)
{
	int i, j, single;

	for (i = 0; i < V; i++) //for each vertex
	{
		single = 0;
		for (j = 0; j < V; j++) // all the vertices
		{
			single += graph[j][i]; //pass all over the the ver in the graph and check if there is a edge with the i vertex (sum the 0 / 1 in the mat, 0=no edge 1=else)
			if (single == 1) //there is edge from the i vertex
			{
				break;
			}
			if ((j == (V - 1)) && (single == 0)) //check all and the sum is 0 (no edges)
			{
				return 1; // is isolated
			}
		}
	}
	return 0; // no exist isolated vertex
}

int connectivity(int** graph)
{
	int i;

	BFS(graph, 1); //call bfs with the first vertex
	for (i = 1; i < V; i++) //check all the others vertices
	{
		if (distance[i] == -1) //check if there is a ver with distance = infinity => unconnected graph
		{
			return 0; //unconnected
		}
	}
	return 1; //connected
}

//**Auxiliary function**//

void BFS(int ** graph, int v) //v=start ver
{
	int i, j;
	dist = 0;
	front = -1, rear = -1;

	for (j = 0; j < V; j++) // reset lists 
	{
		queue[j] = 0; //no ver in the Q
		visited[j] = 0; //no visited
		parent[j] = -1; //NULL
		distance[j] = -1; //infinity
	}


	push_queue(v); // add ver first to q

	while (!isEmpty_queue()) // while the q is not empty
	{
		v = pop_queue(); //pick the first vertex
		visited[v - 1] = 1;

		for (i = 0; i < V; i++)
		{
			if (visited[i])   //if it has already been visited by some other neighbors vertex, it should not be printed again.
			{
				continue;
			}
			if ((graph[v - 1][i] == 1) && (visited[i] == 0)) //If there is 1 and no visited
			{
				push_queue(i + 1); //add vertex to Q
				parent[i] = v; //update parent of each ver
				visited[i] = 1; //update that we visit this vertex
				if (distance[v - 1] == -1) // if its the "root" so we cant do deg +1 so we update the distanse=1
				{
					distance[i] = 1;
				}
				else
				{
					distance[i] = distance[v - 1] + 1; //update distance to be the parent ver distance + 1
				}
				if (dist < distance[i]) //update the distance from the S vertex
				{
					dist = distance[i];
				}
			}
		}
	}
}

void push_queue(int vertex) //ADD VERTEX TO QUEUE
{
	if (rear == V + 1)
		printf("Queue Overflow\n");
	else
	{
		if (front == -1)
		{
			front = 0; //front = next
		}
		rear = rear + 1; //reat = next
		queue[rear] = vertex; //add ver to Q
	}
}

int isEmpty_queue() //CHECK IF THE Q IS EMPTY
{
	if (front == -1 || front > rear) // if front is in the start (no ver) or if front bigger the the rear (pass the end)
	{
		return 1; //empty Q
	}
	else
	{
		return 0; //no empty Q
	}
}

int pop_queue() // DELETE FROM THE FRONT OF Q
{
	int delete_item;
	if (front == -1 || front > rear)
	{
		printf("Queue Underflow\n");
		exit(1);
	}

	delete_item = queue[front]; // get the first ver in the Q
	front = front + 1; //front = next
	return delete_item; //return the first ver in the Q
}

void freeDynamicMatrix(int **mat, int n) // Dynamic 2D array free function
{
	int i;
	for (i = 0; i < n; i++)
	{
		free(mat[i]); // free each line
	}
	free(mat);
}





