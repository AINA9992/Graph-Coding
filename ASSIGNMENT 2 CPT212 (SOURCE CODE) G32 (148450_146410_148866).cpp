/*
AQILAH SYAZANI BINTI NOORAZIZ   _   146410
AINA SYARAFINA BINTI ROSLAN     _   148450
ZAKIRAH BINTI SHARIFUDDIN       _   148866
*/

#include <windows.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <bits/stdc++.h>
#include <limits.h>
using namespace std;

void menu1();//main
void menu2();//sub menu
void menu3();//menu shortest
void exit();
int distance(int, int);
 
// A utility function to add an edge in a directed graph.
void addEdge(vector<int> adj[], int u, int v)
{
    adj[u].push_back(v);
}
 
// A utility function to delete an edge in a directed graph.
void delEdge(vector<int> adj[], int u, int v)
{
    // Traversing through the first vector list and removing the second element from it
    for (int i = 0; i < adj[u].size(); i++) {
        if (adj[u][i] == v) {
            adj[u].erase(adj[u].begin() + i);
            break;
        }
    }
}
 
// A utility function to print the adjacency list representation of graph
void printGraph(vector<int> adj[], int V)
{
	cout << "\n\n\n\t\t ============================\n";
	cout << "\t\t\tDEFAULT GRAPH";
	cout << "\n\t\t ============================\n";
    for (int v = 0; v < V; ++v) {
        cout << "\t\t     vertex " << v << " ";
        for (auto x : adj[v])
            cout << "-> " << x;
        printf("\n");
    }
    printf("\n");
}
 
// A utility function to check and add an edge in a directed graph.
void check_addEdge(vector<int> adj[], int u, int v)
{
	if (adj[u].size() != 0)
	{
		int q = 1;
		// Traversing through the first vector list and add the second element from it
	    for (int i = 0; i < adj[u].size(); i++) 
		{
	        if (adj[u][i] == v) 
			{
				cout << "\n\n\n\tEdge is already available.";
				q = 2;
	            break;
	        }
	    }
	    
	    if (q == 1) 
		{
			cout << "\n\n\n\tEdge created.";
		  	addEdge(adj, u, v); 
	    }
	}
	
	//for vertices that did not have any edges
	else
	{
		cout << "\n\n\n\tEdge created.";
		addEdge(adj, u, v); 
	}	  
	
}

class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // An array of adjacency lists
 	
    // A recursive function to print DFS starting from v
    void BFSUtil(int v, bool visited[]);
    bool isCyclicUtil(int m, bool visited[], bool *rs);  // used by isCyclic()
    
public:
 
    // Constructor and Destructor
    Graph(int V) { this->V = V;  adj = new list<int>[V];}
    ~Graph() { delete [] adj; }
 
    // Method to add an edge
    void addEdge(int v, int w);
 
    // The main function that returns true if the
    // graph is strongly connected, otherwise false
    bool isSC();
 
    // Function that returns reverse (or transpose)
    // of this graph
    Graph getTranspose();
    bool isCyclic();    // returns true if there is a cycle in this graph
};
 
// A recursive function to print DFS starting from v
void Graph::BFSUtil(int v, bool visited[])
{
    // Create a queue for BFS
    list<int> queue;
 
    // Mark the current node as visited and enqueue it
    visited[v] = true;
    queue.push_back(v);
 
    // 'i' will be used to get all adjacent vertices
    // of a vertex
    list<int>::iterator i;
 
    while (!queue.empty())
    {
        // Dequeue a vertex from queue
        v = queue.front();
        queue.pop_front();
 
        // Get all adjacent vertices of the dequeued vertex s
        // If a adjacent has not been visited, then mark it
        // visited and enqueue it
        for (i = adj[v].begin(); i != adj[v].end(); ++i)
        {
            if (!visited[*i])
            {
                visited[*i] = true;
                queue.push_back(*i);
            }
        }
    }
}
 
// Function that returns reverse (or transpose) of this graph
Graph Graph::getTranspose()
{
    Graph g(V);
    for (int v = 0; v < V; v++)
    {
        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;
        for (i = adj[v].begin(); i != adj[v].end(); ++i)
            g.adj[*i].push_back(v);
    }
    return g;
}

// This function is a variation of DFSUtil() in https://www.geeksforgeeks.org/archives/18212
bool Graph::isCyclicUtil(int v, bool visited[], bool *recStack)
{
 	if(visited[v] == false)
    {
        // Mark the current node as visited and part of recursion stack
        visited[v] = true;
        recStack[v] = true;
 
        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;
        for(i = adj[v].begin(); i != adj[v].end(); ++i)
        {
            if ( !visited[*i] && isCyclicUtil(*i, visited, recStack) )
                return true;
            else if (recStack[*i])
                return true;
        }
 
    }
    recStack[v] = false;  // remove the vertex from recursion stack
    return false;
}
 
// Returns true if the graph contains a cycle, else false.
// This function is a variation of DFS() in https://www.geeksforgeeks.org/archives/18212
bool Graph::isCyclic()
{
    // Mark all the vertices as not visited and not part of recursion
    // stack
    bool *visited = new bool[V];
    bool *recStack = new bool[V];
    for(int i = 0; i < V; i++)
    {
        visited[i] = false;
        recStack[i] = false;
    }
 
    // Call the recursive helper function to detect cycle in different
    // DFS trees
    for(int i = 0; i < V; i++)
        if (isCyclicUtil(i, visited, recStack))
            return true;
 
    return false;
} 
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
}
 
// The main function that returns true if graph
// is strongly connected
bool Graph::isSC()
{
    // St1p 1: Mark all the vertices as not
    // visited (For first BFS)
    bool visited[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
 
    // Step 2: Do BFS traversal starting
    // from first vertex.
    BFSUtil(0, visited);
 
    // If BFS traversal doesn’t visit all
    // vertices, then return false.
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
             return false;
 
    // Step 3: Create a reversed graph
    Graph gr = getTranspose();
 
    // Step 4: Mark all the vertices as not
    // visited (For second BFS)
    for(int i = 0; i < V; i++)
        visited[i] = false;
 
    // Step 5: Do BFS for reversed graph
    // starting from first vertex.
    // Staring Vertex must be same starting
    // point of first DFS
    gr.BFSUtil(0, visited);
 
    // If all vertices are not visited in
    // second DFS, then return false
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
             return false;
 
    return true;
}

struct AdjListNode
{
    int dest;
    int weight;
    struct AdjListNode* next;
};

// A structure to represent
// an adjacency list
struct AdjList
{
     
   // Pointer to head node of list
   struct AdjListNode *head;
};
 
// A structure to represent a graph.
// A graph is an array of adjacency lists.
// Size of array will be V (number of
// vertices in graph)
struct Graph5
{
    int V;
    struct AdjList* array;
};

// A utility function to create
// a new adjacency list node
struct AdjListNode* newAdjListNode(
                   int dest, int weight)
{
    struct AdjListNode* newNode =
            (struct AdjListNode*)
      malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

// A utility function that creates
// a graph of V vertices
struct Graph5* createGraph(int V)
{
    struct Graph5* graph5 = (struct Graph5*)
            malloc(sizeof(struct Graph5));
    graph5->V = V;
 
    // Create an array of adjacency lists. 
    // Size of array will be V
    graph5->array = (struct AdjList*)
       malloc(V * sizeof(struct AdjList));
 
    // Initialize each adjacency list
    // as empty by making head as NULL
    for (int i = 0; i < V; ++i)
        graph5->array[i].head = NULL;
 
    return graph5;
}

// Adds an edge to an undirected graph
void addEdge(struct Graph5* graph5, int src,
                   int dest, int weight)
{
    // Add an edge from src to dest. 
    // A new node is added to the adjacency
    // list of src.  The node is
    // added at the beginning
    struct AdjListNode* newNode =
            newAdjListNode(dest, weight);
    newNode->next = graph5->array[src].head;
    graph5->array[src].head = newNode;
 
    // Since graph is undirected,
    // add an edge from dest to src also
//    newNode = newAdjListNode(src, weight);
//    newNode->next = graph5->array[dest].head;
//    graph5->array[dest].head = newNode;
}

// Structure to represent a min heap node
struct MinHeapNode
{
    int  v;
    int dist;
};

// Structure to represent a min heap
struct MinHeap
{
     
    // Number of heap nodes present currently
    int size;    
   
    // Capacity of min heap
    int capacity; 
   
    // This is needed for decreaseKey()
    int *pos;   
    struct MinHeapNode **array;
};

// A utility function to create a
// new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v,
                                 int dist)
{
    struct MinHeapNode* minHeapNode =
           (struct MinHeapNode*)
      malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}
 
// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity)
{
    struct MinHeap* minHeap =
         (struct MinHeap*)
      malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(
            capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array =
         (struct MinHeapNode**)
                 malloc(capacity *
       sizeof(struct MinHeapNode*));
    return minHeap;
}

// A utility function to swap two
// nodes of min heap.
// Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a,
                     struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to
// heapify at given idx
// This function also updates
// position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap,
                                  int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;
 
    if (left < minHeap->size &&
        minHeap->array[left]->dist <
         minHeap->array[smallest]->dist )
      smallest = left;
 
    if (right < minHeap->size &&
        minHeap->array[right]->dist <
         minHeap->array[smallest]->dist )
      smallest = right;
 
    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode =
             minHeap->array[smallest];
        MinHeapNode *idxNode =
                 minHeap->array[idx];
 
        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;
 
        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest],
                         &minHeap->array[idx]);
 
        minHeapify(minHeap, smallest);
    }
}
 
// A utility function to check if
// the given minHeap is ampty or not
int isEmpty(struct MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// Standard function to extract
// minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap*
                                   minHeap)
{
    if (isEmpty(minHeap))
        return NULL;
 
    // Store the root node
    struct MinHeapNode* root =
                   minHeap->array[0];
 
    // Replace root node with last node
    struct MinHeapNode* lastNode =
         minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;
 
    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;
 
    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);
 
    return root;
}
 
// Function to decreasy dist value
// of a given vertex v. This function
// uses pos[] of min heap to get the
// current index of node in min heap
void decreaseKey(struct MinHeap* minHeap,
                         int v, int dist)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];
 
    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;
 
    // Travel up while the complete
    // tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist <
           minHeap->array[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] =
                                      (i-1)/2;
        minHeap->pos[minHeap->array[
                             (i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], 
                 &minHeap->array[(i - 1) / 2]);
 
        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v)
{
   if (minHeap->pos[v] < minHeap->size)
     return true;
   return false;
}

// Function to print shortest
// path from source to j
// using parent array
void printPath(int parent[], int j) // j = destination
{
       
    // Base Case : If j is source
    if (parent[j] == - 1)
        return;
   
    printPath(parent, parent[j]);
   
   	cout << "   " << j;
}

int printSolution(int dist[], int n, 
                      int parent[], int strt, int des)
{
	if (dist[des] == 2147483647) 
	{
		cout << "\n\n\n\tPath did not found. Please add new edges! \n ";
	}
	else
	{
		cout << "\n\n\t\t==========================================================";
		cout << "\n\t\t|    Vertex    |     Distance     |         Path         |";
		cout << "\n\t\t==========================================================";
        cout <<"\n\t\t     " << strt << " -> " << des << "            " << dist[des] << "           " << strt ;
        printPath(parent, des);
	}
}


// The main function that calculates
// distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
void dijkstra(struct Graph5* graph5, int src, int des)
{     
    // Get the number of vertices in graph
    int V = graph5->V;
   
    // dist values used to pick
    // minimum weight edge in cut
    int dist[V];    
 	bool sptSet[V];
 	int parent[V];
 	 for (int i = 0; i < V; i++)
    {
        parent[0] = -1;
        dist[i] = INT_MAX;
        sptSet[i] = false;
    }
    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);
 
    // Initialize min heap with all
    // vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v,
                                      dist[v]);
        minHeap->pos[v] = v;
    }
// 	dist[src] = 0;
 	
    // Make dist value of src vertex
    // as 0 so that it is extracted first
    minHeap->array[src] =
          newMinHeapNode(src, dist[src]);
    minHeap->pos[src]   = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);
 
    // Initially size of min heap is equal to V
    minHeap->size = V;
 
    // In the followin loop,
    // min heap contains all nodes
    // whose shortest distance
    // is not yet finalized.
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with
        // minimum distance value
        struct MinHeapNode* minHeapNode =
                     extractMin(minHeap);
       
        // Store the extracted vertex number
        int u = minHeapNode->v;
        
 		sptSet[u] = true;
 		
        // Traverse through all adjacent
        // vertices of u (the extracted
        // vertex) and update
        // their distance values
        struct AdjListNode* pCrawl =
                     graph5->array[u].head;
        
        while (pCrawl != NULL)
        {
            int v = pCrawl->dest;
 
            // If shortest distance to v is
            // not finalized yet, and distance to v
            // through u is less than its
            // previously calculated distance
            if (isInMinHeap(minHeap, v) &&
                      dist[u] != INT_MAX &&
              pCrawl->weight + dist[u] < dist[v])
            {
                dist[v] = dist[u] + pCrawl->weight;
 				
                // update distance
                // value in min heap also
                decreaseKey(minHeap, v, dist[v]);
				
				if (parent[v] != parent [v-1])
				{
					parent[v] = u;
				}
				
            }
			
            pCrawl = pCrawl->next;
            
        }
    }
		printSolution(dist, V, parent, src, des);
}
 
 
 

 
// Driver code
int main()
{
    int V = 5;
    int dist = 0;
    int choice1, choice2, choice3, vertices1, vertices2, random1, random2;
    
    srand(time(0)); 
    
    Graph g1(5);
    vector<int> adj[V];
 
    // Adding edge as shown in the example figure
    addEdge(adj, 0, 1);
    addEdge(adj, 0, 4);
    addEdge(adj, 2, 1);
    addEdge(adj, 3, 2);
    addEdge(adj, 4, 3);
 	
 	struct Graph5* graph5 = createGraph(V);

	do{
		
	menu1(); // function call to display menu
	
    cout<<"\t\tEnter choice : ";
    cin>>choice1;
    				
    while(cin.fail()) // input must same data type as what we had declared //VALIDATE choice
	{ 
        cout<<"\n\n\t\tInvalid input!!!\n\n\t\t";
        cin.clear();
        cin.ignore(256,'\n');
        system("pause");
		system("cls");
		menu1();
		cout<<"\t\tEnter choice : ";
        cin>>choice1;
    }

	cout << "\n\n\n";

	switch(choice1){
		
				case 1: // strong
						system("cls");	
						cout << "\n\n\n";
						
						do{
		
							cout << "\n\n\t\t\t    ==============================\n";
							cout << "                               STRONGLY CONNECTED GRAPH";
							cout << "\n\t\t\t    ==============================\n\n";
							cout << "                        | 1 | Is the Graph Strongly Connected?|\n";
							
							menu2(); // function call to display menu
							
						    cout<<"\t\tEnter choice : ";
						    cin>>choice2;
						    				
						    while(cin.fail()) // input must same data type as what we had declared //VALIDATE choice
							{ 
						        cout<<"\n\n\t\tInvalid input!!!\n\n\t\t";
						        cin.clear();
						        cin.ignore(256,'\n');
						        system("pause");
								system("cls");
								cout << "\n\n\t\t\t    ==============================\n";
								cout << "                               STRONGLY CONNECTED GRAPH";
								cout << "\n\t\t\t    ==============================\n\n";
								cout << "                        | 1 | Is the Graph Strongly Connected?|\n";
								menu2();
								cout<<"\t\tEnter choice : ";
						        cin>>choice2;
						    }
						
							cout << "\n\n\n";
						
							switch(choice2){
								
										case 1: // check
												system("cls");	
												cout << "\n\n\n\tIs the graph is strongly connected?  :  ";
												
												//copy masuk list
												for (int v = 0; v < V; ++v) 
												{
											        for (auto x : adj[v])
											        {
											        	g1.addEdge(v, x);
													}
												}
												
												g1.isSC()? cout << "Yes\n" : cout << "No\n";
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
															
						 				case 2: // display all edge
												system("cls");	
												cout << "\n\n\n";
												printGraph(adj, V);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
						 						
											    printf("\n");
						 						
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 3: // add random edge
												system("cls");	
												cout << "\n\n\n";
												random1 = (rand() % 5 ) + 1; 
												random2 = (rand() % 5) + 1; 
												random1 = random1 - 1;
												random2 = random2 - 1;
//												
												if (random1 != random2)
												{
													cout << "\n\n\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2 << endl;
													check_addEdge(adj, random1, random2);
													
							 						cout << "\n\n\n\tCompleted..!\n\n\t";
												}
												
												else if (random1 == random2)
												{
													cout << "\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2 << endl;
													cout << "\n\n\n\tCannot create edge..!\n\n\t";
												}
												
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 4: //remove an edge
						 						{	
												system("cls");	
												cout << "\n\n\n";
												cout << "\n\n\n\tEnter vertices 1 : ";
												cin >> vertices1;
												cout << "\n\tEnter vertices 2 : ";
												cin >> vertices2;
												delEdge(adj, vertices1, vertices2);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
												}
						 						
						 				case 5: // back to main menu
						 						break;
						 				default: 
										 		cout<<"\n\n\t\tChoice not valid...\n\n\t\t";
										 		system("pause");
												system("cls");
									 }
							}while (choice2 != 5);
						
 						cout << "\n\n\n\tCompleted..!\n\n\t";
//						system("pause");
						system("cls");
 						break;
									
 				case 2: // cycle
						system("cls");	
						cout << "\n\n\n";
						
						do{
		
							cout << "\n\n\t\t\t\t====================\n";
							cout << "                                    CYCLIC GRAPH";
							cout << "\n\t\t\t\t====================\n\n";
							cout << "                        | 1 | Is the Graph have Cycle?        |\n";
							
							menu2(); // function call to display menu
							
						    cout<<"\t\tEnter choice : ";
						    cin>>choice2;
						    				
						    while(cin.fail()) // input must same data type as what we had declared //VALIDATE choice
							{ 
						        cout<<"\n\n\t\tInvalid input!!!\n\n\t\t";
						        cin.clear();
						        cin.ignore(256,'\n');
						        system("pause");
								system("cls");
								cout << "\n\n\t\t\t\t====================\n";
								cout << "                                    CYCLIC GRAPH";
								cout << "\n\t\t\t\t====================\n\n";
								cout << "                        | 1 | Is the Graph have Cycle?        |\n";
								menu2();
								cout<<"\t\tEnter choice : ";
						        cin>>choice2;
						    }
						
							cout << "\n\n\n";
						
							switch(choice2){
								
										case 1: // check
												system("cls");	
												cout << "\n\n\n";
												for (int v = 0; v < V; ++v) 
												{
											        for (auto x : adj[v])
											        {
											        	g1.addEdge(v, x);
													}
												}
												if(g1.isCyclic())
										        cout << "\n\n\n\tGraph contains cycle";
										    	else
										        cout << "\n\n\n\tGraph doesn't contain cycle";
										
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
															
						 				case 2: // display all edge
												system("cls");	
												cout << "\n\n\n";
												printGraph(adj, V);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 3: // add random edge
												system("cls");	
												cout << "\n\n\n";

												random1 = (rand() % 5 ) + 1; 
												random2 = (rand() % 5) + 1; 
												random1 = random1 - 1;
												random2 = random2 - 1;
												
												if (random1 != random2)
												{
													cout << "\n\n\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2 << endl;
													check_addEdge(adj, random1, random2);
													
							 						cout << "\n\n\n\tCompleted..!\n\n\t";
												}
												
												else if (random1 == random2)
												{
													cout << "\n\n\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2;
													cout << "\n\n\n\tCannot create edge..!\n\n\t";
												}
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 4: //remove an edge
						 						{	
												system("cls");	
												cout << "\n\n\n";
												cout << "\n\n\n\tEnter vertices 1 : ";
												cin >> vertices1;
												cout << "\n\tEnter vertices 2 : ";
												cin >> vertices2;
												delEdge(adj, vertices1, vertices2);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
												}
						 						
						 				case 5: // back to main menu
						 						break;
						 				default: 
										 		cout<<"\n\n\t\tChoice not valid...\n\n\t\t";
										 		system("pause");
												system("cls");
									 }
							}while (choice2 != 5);
						
 						cout << "\n\n\n\tCompleted..!\n\n\t";
						system("cls");
 						break;
 						
 				case 3: // shortest path
						system("cls");	
						cout << "\n\n\n";
						
						do{
								
							menu3(); // function call to display menu
							
						    cout<<"\t\tEnter choice : ";
						    cin>>choice3;
						    				
						    while(cin.fail()) // input must same data type as what we had declared //VALIDATE choice
							{ 
						        cout<<"\n\n\n\t\tInvalid input!!!\n\n\t\t";
						        cin.clear();
						        cin.ignore(256,'\n');
						        system("pause");
								system("cls");
								menu3();
								cout<<"\t\tEnter choice : ";
						        cin>>choice3;
						    }
						
							cout << "\n\n\n";
						
							switch(choice3){
								
										case 1: // insert two vertices
												system("cls");	
												cout << "\n\n\n";
												cout << "\n\n\n\tEnter vertices 1 : ";
												cin >> vertices1;
												cout << "\n\tEnter vertices 2 : ";
												cin >> vertices2;
												
												for (int v = 0; v < V; ++v) 
												{
											        for (auto x : adj[v])
											        {
											        	dist = distance(v, x);
											        	addEdge(graph5, v, x, dist);
													}
												}
												
												// calculate
												dijkstra(graph5, vertices1, vertices2);
												
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 2: // display all edge
												system("cls");	
												cout << "\n\n\n";
												printGraph(adj, V);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 3: // add random edge
												system("cls");	
												cout << "\n\n\n";

												random1 = (rand() % 5 ) + 1; 
												random2 = (rand() % 5) + 1; 
												random1 = random1 - 1;
												random2 = random2 - 1;
												
												if (random1 != random2)
												{
													cout << "\n\n\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2 << endl;
													check_addEdge(adj, random1, random2);
													
							 						cout << "\n\n\n\tCompleted..!\n\n\t";
												}
												
												else if (random1 == random2)
												{
													cout << "\n\n\n\tRandom vertices 1 : " << random1;
													cout << "\n\tRandom vertices 2 : " << random2;
													cout << "\n\n\n\tCannot create edge..!\n\n\t";
												}
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
						 						
						 				case 4: //remove an edge
						 						{	
												system("cls");	
												cout << "\n\n\n";
												cout << "\n\n\n\tEnter vertices 1 : ";
												cin >> vertices1;
												cout << "\n\tEnter vertices 2 : ";
												cin >> vertices2;
												delEdge(adj, vertices1, vertices2);
						 						cout << "\n\n\n\tCompleted..!\n\n\t";
												system("pause");
												system("cls");
						 						break;
												}
						 						
						 				case 5: // back to main menu
						 						break;
						 				default: 
										 		cout<<"\n\n\t\tChoice not valid...\n\n\t\t";
										 		system("pause");
												system("cls");
									 }
							}while (choice3 != 5);
						
						system("cls");
 						break;
 						
 				case 4: //reset
 						system("cls");	
						cout << "\n\n\n";
						
						for (int i = 0; i < 5; i++)
						{
							for (int j = 0; j < 5; j++)
							{
								delEdge(adj, i, j);
							}
						}
						
						addEdge(adj, 0, 1);
					    addEdge(adj, 0, 4);
					    addEdge(adj, 2, 1);
					    addEdge(adj, 3, 2);
					    addEdge(adj, 4, 3);
					    
					    g1.addEdge(0, 1);
					    g1.addEdge(0, 4);
					    g1.addEdge(2, 1);
					    g1.addEdge(3, 2);
					    g1.addEdge(4, 3);
					    
 						cout << "\n\n\n\tReset Completed..!\n\n\t";
						system("pause");
						system("cls");
 						break;
 						
 				case 5: //end session
		 				exit();
 						break;
 				default: 
				 		cout<<"\n\n\t\tChoice not valid...\n\n\t\t";
				 		system("pause");
						system("cls");
			 }
	}while (choice1 != 5);
 
    return 0;
}

void menu1()
{
	cout << "\n\n\t\t\t\t====================\n";
	cout << "                                        MENU";
	cout << "\n\t\t\t\t====================\n\n";
	cout << "                        | 1 | Strongly Connected Graph    |\n";
	cout << "                        | 2 | Cyclic Graph                |\n";
	cout << "                        | 3 | Find the Shortest Path      |\n";
	cout << "                        | 4 | Reset Graph                 |\n";
	cout << "                        | 5 | End Session                 |\n\n\n\n\n\t";
}

void menu2()
{
	cout << "                        | 2 | Display Available Edges         |\n";
	cout << "                        | 3 | Add Random Edge                 |\n";
	cout << "                        | 4 | Remove an Edge                  |\n";
	cout << "                        | 5 | Back to the Main Menu           |\n\n\n\n\n\t";
}

void menu3()
{
	cout << "\n\n\t\t\t\t====================\n";
	cout << "                                    SHORTEST PATH";
	cout << "\n\t\t\t\t====================\n\n";
	cout << "                        | 1 | Compute Shortest Path?          |\n";
	cout << "                        | 2 | Display Available Edges         |\n";
	cout << "                        | 3 | Add Random Edge                 |\n";
	cout << "                        | 4 | Remove an Edge                  |\n";
	cout << "                        | 5 | Back to the Main Menu           |\n\n\n\n\n\t";
}

void exit() // display exit note
{
	system("cls");
	cout << "\n\n\n\n\n\n\n\t\t**  *     *    *    *   *   *   *       *    *   *    *     **\n";
	cout << "\t\t   *     *     *   *     *   **  *   *  *         *  *   *     *   *     *\n";
	cout << "\t\t   *     *   *   * * *   *             *     *     *   *     *\n";
	cout << "\t\t   *     *     *   *     *   *  *   *  *           *     *     *   *     *\n";
	cout << "\t\t   *     *     *   *     *   *   *   *   *          *      *     **\n";
	Sleep (50);
	cout<< "\n\n\t______\n\n\n";
	cout<< "\t\t\t\t\t\t  for using this program\n\n";
	cout << "\t\t\t\t\t       AND PLEASE COME AGAIN...  :)";
	cout<< "\n\n\t______\n\n\n";
}

int distance(int x, int y)
{
	if (x == 0)
	{
		if (y == 1)
			return 7672;
		else if (y == 2)
			return 8381;
		else if (y == 3)
			return 13571;
		else if (y == 4)
			return 9737;
	}
	else if ( x == 1)
	{
		if (y == 0)
			return 7672;
		else if (y == 2)
			return 915;
		else if (y == 3)
			return 7466;
		else if (y == 4)
			return 7962;
	}
	else if ( x == 2)
	{
		if (y == 0)
			return 8381;
		else if (y == 1)
			return 915;
		else if (y == 3)
			return 7466;
		else if (y == 4)
			return 7962;
	}
	else if ( x == 3)
	{
		if (y == 0)
			return 13751;
		else if (y == 1)
			return 7466;
		else if (y == 2)
			return 7317;
		else if (y == 4)
			return 4853;
	}
	else if ( x == 4)
	{
		if (y == 0)
			return 9737;
		else if (y == 1)
			return 7962;
		else if (y == 2)
			return 8506;
		else if (y == 3)
			return 4853;
	}
}
