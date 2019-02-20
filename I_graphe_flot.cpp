#include "I_graphe_flot.h"
#include <sys/time.h>
#include <queue>


void Barf(const char *s)
{
  cerr << s << endl;
  exit(-1);
}


char *Alloc(int n)
{
  char *p;
  p = new char[n];

  if (p == NULL)
    Barf("Out of space");

  return p;
}


double Min(double x,double y)
{
  return (x > y) ? y : x;
}

double Max(double x,double y)
{
  return (x > y) ? x : y;
}

double Abs(double x)
{
  return (x > 0) ? x : -x;
}

void InitRandom (int seed)
{
    struct timeval tp;

    if (seed == 0){
        gettimeofday(&tp, 0);
        srandom(tp.tv_sec + tp.tv_usec);
    }
    else
	srandom(seed);
}

Queue *MakeQueue(int n)
{
  Queue *Q;

  Q = (Queue *) Alloc(sizeof(Queue));

  Q->data = (int *) Alloc(n * sizeof(int));
  Q->tail = 0;
  Q->head = 0;

  Q->size = n;

  return Q;
}


int Dequeue(Queue *Q)
{
  int v;

  if (Q->tail == Q->head)
    Barf("Attempt to dequeue from empty queue");

  v = Q->data[Q->head];
  Q->head = (Q->head == Q->size - 1) ? 0 : Q->head + 1;

  return v;
}

void Enqueue(Queue *Q,int k)
{
  if (Q->head == Q->tail + 1 ||
      (Q->tail == Q->size - 1 && Q->head == 0))
    Barf("Queue overfull");

  Q->data[Q->tail] = k;
  Q->tail = (Q->tail == Q->size - 1) ? 0 : Q->tail + 1;
}

bool QueueEmpty(Queue *Q)
{
  return(Q->head == Q->tail);
}



/************  Classe I_Graph_Flot  ***********************/

/* for file output */
void I_Graph_Flot::OutputFlow(const char *fich, double s)
{
  	int i;
  	ofstream f(fich);

  	f << "s " << s << endl;
	for (i = 0; i <= max_v; i++)
  	{
  		WriteVertex2(i, f);
  	}
  	f.close();
}


void I_Graph_Flot::WriteVertex(int v, ofstream &f)
{
  	I_Edge_Flot *e;

  	e = A[v];

	while (e != (I_Edge_Flot *) 0)
	{
    		if (e->c > 0)
		{
       			f << "a " << e->t << " " << e->h << " " << e->c << endl;
    		}
    		e = e->next;
  	}
}


void I_Graph_Flot::WriteVertex2(int v, ofstream &f)
{
  	I_Edge_Flot *e;

  	e = A[v];

	while (e != (I_Edge_Flot *) 0)
	{
    		if (e->f > 0)
		{
			f << "f " << e->t << " " << e->h << " " << e->f << endl;
    		}
   		 e = e->next;
  	}
}

void I_Graph_Flot::WriteVertex3(int v, ofstream &f)
{
  	I_Edge_Flot *e;

  	e = A[v];

	while (e != (I_Edge_Flot *) 0)
	{
    		if (e->c > 0)
		{
        		f << e->t << " " << e->h << " " << e->c << " " << e->f << endl;
    		}
    		e = e->next;
  	}
}

bool I_Graph_Flot::InputFlowGraph(char *fich,int &s,int &t)
{
  	int i;
  	char buff[256];
  	int nodes, edges, u, v;
  	double cap;

  	ifstream f(fich);

	if(!f)
	{
		cout << "Impossible d'ouvrir le fichier " << fich << endl;
		return false;
	}


  	f >> buff;
  	f >> buff;
  	f >> nodes;
  	f >> edges;

  	f >> buff;
  	f >> s;
  	f >> buff;

  	f >> buff;
  	f >> t;
  	f >> buff;

  	InitGraph(nodes);

  	for(i=0;i<edges;i++)
  	{
   		f >> buff;
    		f >> u;
    		f >> v;
    		f >> cap;

    		AddEdge(u-1,v-1,0,cap);
  	}

  	for (i = 0; i < nodes; i++)
    		AddVertex(i);

	f.close();

	return true;
}


void I_Graph_Flot::InputFlowGraph(list<simpleEdge> &Te,int nodes,int &s,int &t)
{
  	int i;
  	int u, v;
  	double cap_inf,cap_sup;

	InitGraph(nodes);

  	list<simpleEdge>::iterator it = Te.begin();

  	while(it != Te.end())
  	{
    	u = (*it).node1;
    	v = (*it).node2;
    	cap_inf = 0;
   		cap_sup = (*it).cap;

    	AddEdge(u-1,v-1,cap_inf,cap_sup);

    	it++;
  	}

  	for (i = 0; i < nodes; i++)
   		AddVertex(i);
}

void I_Graph_Flot::InitGraph(int n)
{
 	int i;

	A = new I_Edge_Flot *[n];
	V = new int[n];

  	for (i = 0; i < n; i++)
	{
    		A[i] = (I_Edge_Flot *) 0;
    		V[i] = FALSE;
  	}

  	size = 0;
  	max_v = -1;
  	edge_count = 0;
}


void I_Graph_Flot::AddVertex(int v)
{
  	if (V[v] == TRUE)
    		Barf("Vertex already present");

  	V[v] = TRUE;
  	size++;

  	if (v > max_v)
    		max_v = v;
}

void I_Graph_Flot::AddEdge(int v1,int v2,double a_inf,double a_sup,int corresp)
{
  	I_Edge_Flot *e1, *e2;

  	if (v1 == v2)
    		Barf("No Loops");

  	if ((e1 = EdgeLookup(v1, v2)) != (I_Edge_Flot *) 0)
  	{
    		e1->c += a_sup;
    		e1->c_inf += a_inf;
    		return;
  	}

  	e1 = new I_Edge_Flot;
	e2 = new I_Edge_Flot;

  	e1->mate = e2;
  	e2->mate = e1;

  	e1->next = A[v1];
  	A[v1] = e1;
  	e1->t = v1;
  	e1->h = v2;
  	e1->c = a_sup;
  	e1->c_inf = a_inf;
  	e1->f = 0;
	//e1->label = 0;
	e1->numDiarc = corresp;

  	e2->next = A[v2];
  	A[v2] = e2;
  	e2->t = v2;
  	e2->h = v1;
  	e2->c = 0;
  	e2->c_inf = 0;
  	e2->f = 0;
	//e2->label = 0;
	e2->numDiarc = corresp;

  	edge_count++;
}

I_Edge_Flot *I_Graph_Flot::EdgeLookup(int v1,int v2)
{
  	I_Edge_Flot *e;

  	e = A[v1];
  	while (e != (I_Edge_Flot *) 0)
  	{
    		if (e->h == v2)
      		return e;

		e = e->next;
  	}

  	return (I_Edge_Flot *) 0;
}


void I_Graph_Flot::GraphViz(char *fichier)
{
	ofstream fich(fichier);
	I_Edge_Flot *e;

	fich << "digraph G {"<< endl;
	fich << "node [shape=circle];" << endl;

	for(int i=0;i<size;i++)
	{
		fich << i << ";" << endl;
	}


	for(int i=0;i<size;i++)
	{
		e = A[i];
		while(e != ((I_Edge_Flot *) 0))
		{
			if((e->c != 0) /*&& (e->f != 0)*/)
			{
				fich << e->t << " -> " << e->h << "[label=\"" << "(";
				fich << e->f << "," << e->c << ")" << "\"];" << endl;
			}
			e = e->next;
		}
	}

	fich << "}" << endl;

	fich.close();
}


void I_Graph_Flot::GraphViz(int s, int t,char *fichier)
{
	ofstream fich(fichier);
	I_Edge_Flot *e;

	fich << "I_Graph_Flot G {"<< endl;
	fich << "node [shape=circle];" << endl;

	fich << s << " [label=\" s=" << s << "\" " << "color = red];" << endl;
	fich << t << " [label=\" t=" << t << "\" " << "color = blue];" << endl;

	for(int i=0;i<size;i++)
	{
		if((i != s) && (i != t))
		fich << i << ";" << endl;
	}


	for(int i=0;i<size;i++)
	{
		e = A[i];
		while(e != ((I_Edge_Flot *) 0))
		{
			if((e->c != 0) /*&& (e->f != 0)*/)
			{
				fich << e->t << " -> " << e->h << "[label=\"" << "(" << setprecision(2)
				<< e->f << "," << e->c_inf << "," << e->c << ")" << "\"];" << endl;
			}

			e = e->next;
		}
	}

	fich << "}" << endl;

	fich.close();

}

void I_Graph_Flot::DeleteGraph_Flot()
{
	I_Edge_Flot *e;

	for(int i=0;i<size;i++)
	{
		while(A[i] != ((I_Edge_Flot *)0))
		{
			e = A[i];
			A[i] = A[i]->next;
			delete e;
		}
	}

	delete[] A;
	delete[] V;
}


void I_Graph_Flot::CorrectionFlot()
{
	I_Edge_Flot *e;

	for(int i=0;i<size;i++)
	{
		e = A[i];
		while(e != ((I_Edge_Flot *) 0))
		{
			if(e->f < 0)
			{
				e->f = 0;
			}

			e = e->next;
		}
	}
}

double I_Graph_Flot::AugmenteFlot(int s,int t,bool *&marquer,bool avecLabel)
{
	queue<int> Q;
	I_Edge_Flot **pred;	//Tableau de prédécesseurs
	int u;
	I_Edge_Flot *e;
	double augm;

	//marquer = new bool[size];
	pred = new I_Edge_Flot*[size];

	for(int i=0;i<size;i++)
	{
		marquer[i] = false;
		pred[i] = (I_Edge_Flot *)0;
	}

	/*cout << "S == " << s << endl;
	cout << "T == " << t << endl;*/

	/*//Initialisation du nombre de label et du parametre
	nbLabel = 0;
	parLabel = edge_count;*/

	marquer[s] = true;
	Q.push(s);
	pred[s] = (I_Edge_Flot *)0;

	while(!Q.empty() && !marquer[t]) //Tant que la file est non vide et que t n'est pas dans la file
	{
		//On sort le premier element de la file
		u = Q.front();
		Q.pop();

		//cout << "U == " << u << endl;

		//on enfile tout ses voisins v tels que les arcs (u,v) ne sont pas saturés
		//et tous les arcs (v,u) tels que le flot est > à la borne inf
		e = A[u];
		while(e != ((I_Edge_Flot *)0))
		{
			//cout << "v == " << e->h << endl;

			//On verifie si pour l'arc qui sort de u, le sommet v
			//est non marqué et si le flot est < la capacité
			if((!marquer[e->h]) && (e->f < e->c))
			{
				Q.push(e->h);
				marquer[e->h] = true;
				pred[e->h] = e;

				//cout << "Pris " << e->h << endl;
			}

			//cout << "mate v == " << e->mate->h << endl;

			//On regarde aussi si pour l'arc entrant en u, le sommet v
			//est non marqué et le flot est > C_inf

			if((!marquer[e->mate->t]) && (e->mate->f > e->mate->c_inf))
			{
				Q.push(e->mate->t);
				marquer[e->mate->t] = true;
				pred[e->mate->t] = e->mate;

				//cout << "Mate pris " << e->mate->t << endl;
			}

			e = e->next;
		}
	}

	//Si on trouve une chaine augmentante
	if(marquer[t])
	{
		augm = INFINIT;
		u = t;
		do
		{
			//cout << u << " -> ";
			e = pred[u];

			if(e->h == u)	//Dans ce cas, l'arc est pris dans le sens direct
			{		//on augmente le flot
				augm = Min(augm,e->c - e->f);
				u = e->t;
			}
			else		//Dans ce cas, l'arc est pris dans le sens inverse
			{		//on diminue le flot

				augm = Min(augm,e->f - e->c_inf);
				u = e->h;
			}
		}while(u != s);

		//cout << s << endl;

		//cout << "Augmentation = " << augm << endl;

		//On met à jour la valeur du flot sur le cycle ainsi crée
		u = t;
		do
		{
			e = pred[u];
			if(e->h == u)	//Arc direct. On augmente le flot
			{
				e->f += augm;
				u = e->t;

				/*//On labelise les arcs du chemin
				nbLabel++;
				e->label = e->label*parLabel + nbLabel;*/
			}
			else	//Arc inverse. On diminue le flot
			{
				e->f -= augm;
				u = e->h;
			}
		}while(u != s);

		/*delete[] marquer;
		marquer = NULL;*/
	}
	else	//On a pas obtenu de chaine augmentante
	{
		augm = 0;
		//cout << "Pas de chaine" << endl;
	}

	delete[] pred;

	return augm;
}

I_Graph_Flot::~I_Graph_Flot()
{
	DeleteGraph_Flot();
}


/*******    Calcul de flot avec Goldberg *********/

double I_Graph_Flot::FindFlow(int s,int t)
{
 	InitFlow();
  	Goldberg(s, t);
  	EndGoldberg();

  	return VertexFlow(s);
}

double I_Graph_Flot::VertexFlow(int i)
{
  	I_Edge_Flot *e;
  	double flow;

  	flow = 0;
  	e = A[i];
  	while (e != (I_Edge_Flot *) 0)
	{
    		flow += e->f;
    		e = e->next;
  	}

  	return flow;
}

void I_Graph_Flot::ValidFlow(int s,int t)
{
  	int i;
  	I_Edge_Flot *e;

  	for (i = 0; i < size; i++)
  	{
    		e = A[i];
    		while (e != (I_Edge_Flot *) 0)
		{
      			if (e->f != -e->mate->f)
			Barf("Antisymmetry violated");
      			if (e->f > e->c)
			Barf("Capacity violated");
      			e = e->next;
    		}
  	}

  	for (i = 0; i< size; i++)
	{
    		if (i == s || i == t)
    	  	continue;

    		if (VertexFlow(i) != 0)
    		{
      			Barf("Conservation violated");
    		}
  	}

  	if (VertexFlow(s) != -VertexFlow(t))
    		Barf("Network leaks!");
}

void I_Graph_Flot::MarkCut(int u,int C[],int *n)
{
  	int S[MAX_N], h, t, i, v, count;
  	I_Edge_Flot *e;

  	for (i = 0; i < size; i++)
    	C[i] = 0;

  	count = 1;
  	h = t = 0;
  	S[0] = u;
  	C[u] = u;
  	while (t <= h)
	{
    		v = S[t++];
    		e = A[v];
    		while (e != (I_Edge_Flot *) 0)
		{
      			if (C[e->h] == 0 && e->f < e->c)
			{
				count++;
				C[e->h] = 1;
				S[++h] = e->h;
      			}
      			e = e->next;
    		}
  	}
  	*n = count;
}

void I_Graph_Flot::PrintCut(int s)
{
  	int C[MAX_N], n, i, count;
  	I_Edge_Flot *e;

  	MarkCut(s, C, &n);
  	printf("Reachable from source %d\n", n);

  	printf("Cut vertices:\n");

  	count = 0;
  	for (i = 0; i < size; i++)
	{
    		if (C[i] == 0)
      		continue;

    		e = A[i];
    		while (e != (I_Edge_Flot *) 0)
		{
      			if (C[e->h] == 0)
			{
				cout << i << " " << endl;
				if (++count % 10 == 0)
	  			cout << endl;
				break;
      			}
      			e = e->next;
    		}
  	}
  	cout << endl;
}

void I_Graph_Flot::InitFlow()
{
  	int i;
  	I_Edge_Flot *e;

  	for (i = 0; i < size; i++)
	{
    		e = A[i];
    		while (e != (I_Edge_Flot *) 0)
		{
      			e->f = 0;
      			e = e->next;
    		}
  	}
}


void I_Graph_Flot::Goldberg(int s,int t)
{
  	InitGoldberg(s);
  	InitRandom(0);

  	RCount = SCount = UCount = 0;
  	Goldberg1(s, t);

  	cout << "G->SCount " << SCount << " G->UCount " << UCount << " G->RCount " << RCount << endl;
}


void I_Graph_Flot::Goldberg1(int s,int t)
{
  	I_Edge_Flot *e;
  	Queue *Q;
  	int v;
  	int count, checkpoint, LabelCount;

  	SetLabels(t, s);

  	Q = MakeQueue(size);
  	e = A[s];
  	while (e != (I_Edge_Flot *) 0)
	{
    		if (e->c > 0 && e->h != t)
      		{
			Enqueue(Q, e->h);
    			cout << "H = " << e->h << endl;
		}

		e = e->next;
  	}

	cout << "GGGGG111111111" << endl;

  	count = 0;   checkpoint = 0;   LabelCount = 0;

  	while (! QueueEmpty(Q))
	{
		count++;
    		checkpoint++;

    		if (count >= edge_count/2)
		{
      			count = 0;
      			LabelCount++;
      			SetLabels(t, s);
    		}

		/*cout << "Excess == " << endl;

		for(int i=0;i<size;i++)
		{
			cout << "E " << i << " " << Excess[i] << endl;
		}

		cout << endl << "Dist == " << endl;

		for(int i=0;i<size;i++)
		{
			cout << "D " << i << " " << Dist[i] << endl;
		}*/



    		v = Dequeue(Q);

		cout << "v = " << v << endl;

		//cout << "V = " << v << " E(" << v << ") = " << Excess[v] << endl;

		Discharge(v, Q, s, t);
    		if (Excess[v] > 0 && Dist[v] < INFINIT)
      		{
			cout << " --> V = " << v << endl;

			Enqueue(Q, v);
  		}

	}

  	cout << "Relabelings :" << LabelCount << endl;
}


void I_Graph_Flot::InitGoldberg(int s)
{
  	int i;
  	I_Edge_Flot *e;

  	/*Initialisation Excess, Dist, Current*/
  	Excess = new double[size];
  	Dist = new int[size];
  	Current = new I_Edge_Flot*[size];

  	for (i = 0; i < size; i++)
	{
    		Excess[i] = 0;
    		Dist[i] = 0;
    		Current[i] = A[i];
  	}

  	Dist[s] = size;

  	e = A[s];
  	while (e != (I_Edge_Flot *) 0)
	{
    		if (e->c > 0)
		{
      			e->f = e->c;
      			e->mate->f = -e->c;
      			Excess[e->h] += e->c;
    		}
    		e = e->next;
  	}
}

void I_Graph_Flot::SetLabels(int v1,int v2)
{
	//int  S[MAX_N], h, t, i, v;
	int  *S, h, t, i, v;
  	I_Edge_Flot *e, *e1;

	S = new int[size];

  	for (i = 0; i < size; i++)
      		Dist[i] = -1;

  	h = t = 0;
  	S[0] = v1;
  	Dist[v1] = 0;

	while (t <= h)
	{
    		v = S[t++];
    		e = A[v];
    		while (e != (I_Edge_Flot *) 0)
		{
      			e1 = e->mate;
      			if (Dist[e1->t] == -1 && e1->f < e1->c)
			{
				Dist[e1->t] = Dist[v] + 1;
				S[++h] = e1->t;
      			}
      			e = e->next;
    		}
  	}

  	if (Dist[v2] == -1)
	{
    		h = t = 0;
    		S[0] = v2;
    		Dist[v2] = size;
    		while (t <= h)
		{
      			v = S[t++];
      			e = A[v];
      			while (e != (I_Edge_Flot *) 0)
			{
				e1 = e->mate;
				if (Dist[e1->t] == -1 && e1->f < e1->c)
				{
	  				Dist[e1->t] = Dist[v] + 1;
	  				S[++h] = e1->t;
				}
				e = e->next;
      			}
    		}
  	}

	delete[] S;
}



void I_Graph_Flot::Discharge(int v,Queue *Q,int s,int t)
{
  	int d;

  	d = Dist[v];
  	while (Dist[v] == d && Excess[v] > 0)
    		PushRelabel(v, Q, s, t);
}


void I_Graph_Flot::PushRelabel(int v,Queue *Q,int s,int t)
{
  	I_Edge_Flot *e;

  	e = Current[v];
  	if (e->c > e->f  && Dist[v] == Dist[e->h] + 1)
	{
    		cout << "E = " << e->t << " " << e->h << " = " << e->c << " " << e->f << endl;

		if (Excess[e->h] == 0 && e->h != s && e->h != t)
      		{
			cout << "Enqueue = " << e->h << endl;
			Enqueue(Q, e->h);
    		}

		Push(e);
  	}
  	else if (e->next == (I_Edge_Flot *) 0)
	{
    		Current[v] = A[v];
    		Relabel(v);
  	}
  	else
   		Current[v] = e->next;
}

void I_Graph_Flot::Push(I_Edge_Flot *e)
{
  	int v, w;
  	double d;


  	v = e->t;  w = e->h;
  	d = Min(e->c - e->f, Excess[v]);

  	if (d == e->c - e->f)
    		SCount++;
  	else
    		UCount++;


    	cout << "Valeur de Push = " << e->t << " " << e->h << " d = " << d << endl;
    	cout << "e->f = " << e->f << " e->mate->f = " << e->mate->f << endl;

  	e->f += d;
  	e->mate->f -= d;
  	Excess[v] -= d;
  	Excess[w] += d;

	cout << "Excess[" << v << "] = " << Excess[v];
	cout << " " << "Excess[" << w << "] = " << Excess[w] << endl;
}

void I_Graph_Flot::Relabel(int v)
{
  	cout << "Dist " << v << " ";

	int d;
  	double r;
  	I_Edge_Flot *e;

  	d = INFINIT;
  	e = A[v];
  	while (e != (I_Edge_Flot *) 0)
	{
    		r = e->c - e->f;
    		if (r > 0 && Dist[e->h] + 1 < d)
      		d = Dist[e->h] + 1;
    		e = e->next;
  	}
  	Dist[v] = d;
  	RCount++;

	cout << Dist[v] << endl;
}

void I_Graph_Flot::EndGoldberg()
{
	delete[] Excess;
	delete[] Dist;
	delete[] Current;

	Excess = NULL;
	Dist = NULL;
	Current = NULL;
}


//*****************************************
//Calcul de flot avec chaine augmentante

double I_Graph_Flot::flotMaxAugmentation(int s,int t,bool *&marquer)
{
	double augm,flot = 0;

	do
	{
		augm = AugmenteFlot(s,t,marquer);
		flot += augm;

		/*if(augm != 0)
			delete[] marquer;*/

	}while(augm > 0);

	return flot;
}


