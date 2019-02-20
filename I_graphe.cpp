#include "I_graphe.h"
/*#include "fifo.h"
#include "pile.h"*/

//static Bool mark_tab[N_TAB];


bool AllocateGraph(long n,long m,I_Graph *g)
{
	long i;

	g->n_Nodes = n;
	g->m_Edges = m;

	/*Allocation of nodes table*/
	g->Nodes = new I_Node[n];

	if(g->Nodes == NULL)
	{
		cout << "Unable to allocate edges memory" << endl;
		return false;
	}

	/*Allocation of edges table*/
	g->Edges = new I_Edge[2*m];
	if(g->Edges == NULL)
	{
		cout << "Unable to allocate edges memory" << endl;

		delete[] g->Nodes;
		return false;
	}

	for(i=0;i<n;i++)
	{
		(g->Nodes[i]).id = i+1;
		(g->Nodes[i]).first_edge = NIL_E;
	}

	return true;
}


bool DeleteGraph(I_Graph *g)
{
	delete[] g->Nodes;
	delete[] g->Edges;

	return true;
}


bool InitGraph_from_file(char *fich,I_Graph *gr,int *k)
{
	long n,m,i,j;
	char c[2];

	I_Edge *eptr1,*eptr2;
	long node1,node2;
	double cap;
	I_Node *nptr1,*nptr2;

	ifstream f(fich);

	if(!f) //== NULL)
	{
		cout << "Erreur: Impossible d'ouvrir le fichier ";
		cout << fich << endl;

		return false;  /// Instead of defined False
	}

	/*Reading of the order of the problem*/
	f >> c;
	f >> (*k);

	/*Reading nodes and edges numbers*/
	f >> c;

	if(c[0] != 'p')
		return false;

	f >> n;
	f >> m;

	if(!AllocateGraph(n,m,gr))
	{
		return false;
	}


	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0L]);
	eptr2 = &(gr->Edges[m]);

	for(j=0;j<m;j++)
	{
		f >> c;
		f >> node1;
		f >> node2;
		f >> cap;

		cout << c << " " << node1 << " " << node2 << " " << cap << endl;

		if((c[0] != 'e') || (node1 < 1) || (node1 > n) || (node2 < 1) || (node2 > n))
		{
			cout << "Error in the input file";
			DeleteGraph(gr);
			return false;  /// Instead of defined False
		}

		--node1;
		--node2;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j+1;
		eptr2->num = j+1;

		if(nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NIL_E;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if(nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NIL_E;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;
	}

	f.close();

	/*Initialisation des autres champs dans la table des sommets*/
	for(i=0;i<n;i++)
		gr->Nodes[i].terminal = false;


	for(i=0L;i<n-1;i++)
	{
		gr->Nodes[i].next_gr = &(gr->Nodes[i+1]);
		gr->Nodes[i].id_W = 1;
		gr->Nodes[i].next_W = &(gr->Nodes[i]);
		gr->Nodes[i].n_sh = &(gr->Nodes[i]);
		gr->Nodes[i].next_sh = &(gr->Nodes[i]);
	}
	gr->Nodes[n-1].next_gr = &(gr->Nodes[0]);
	gr->Nodes[n-1].id_W = 1;
	gr->Nodes[n-1].next_W = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].n_sh = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].next_sh = &(gr->Nodes[n-1]);


	return true;
}



bool InitGraph_from_list(simpleEdge *Te,long n,long m,I_Graph *gr)
{
	long i,j;
	I_Edge *eptr1,*eptr2;
	long node1,node2;
	double cap;
	I_Node *nptr1,*nptr2;

	if(!AllocateGraph(n,m,gr))
	{
		return false;  /// Instead of defined False
	}

	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0]);
	eptr2 = &(gr->Edges[m]);

	for(j=0L;j<m;j++)
	{
		node1 = Te[j].node1;
		node2 = Te[j].node2;
		cap = Te[j].cap;

		--node1;
		--node2;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j+1;
		eptr2->num = j+1;

		if(nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NULL;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if(nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NULL;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;

	}

	for(i=0;i<n;i++)
		gr->Nodes[i].terminal = false;

	/*Initialisation des autres champs dans la table des sommets*/
	for(i=0L;i<n-1;i++)
	{
		gr->Nodes[i].next_gr = &(gr->Nodes[i+1]);
		gr->Nodes[i].id_W = 1;
		gr->Nodes[i].next_W = &(gr->Nodes[i]);
		gr->Nodes[i].n_sh = &(gr->Nodes[i]);
		gr->Nodes[i].next_sh = &(gr->Nodes[i]);
	}
	gr->Nodes[n-1].next_gr = &(gr->Nodes[0]);
	gr->Nodes[n-1].id_W = 1;
	gr->Nodes[n-1].next_W = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].n_sh = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].next_sh = &(gr->Nodes[n-1]);

	return true;
}


bool InitGraph_from_list(vector<simpleEdge> &Te,long n,I_Graph *gr)
{
	int i,j;
	I_Edge *eptr1,*eptr2;
	long node1,node2;
	double cap;
	I_Node *nptr1,*nptr2;
	int m = Te.size();

	if(!AllocateGraph(n,m,gr))
	{
		return false;  /// Instead of defined False
	}

	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0]);
	eptr2 = &(gr->Edges[m]);


	for(j=0;j<m;j++)
	{
		node1 = Te[j].node1;
		node2 = Te[j].node2;
		cap = Te[j].cap;

		--node1;
		--node2;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j+1;
		eptr2->num = j+1;

		if(nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NULL;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if(nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NULL;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;

	}

	for(i=0;i<n;i++)
		gr->Nodes[i].terminal = false;

	/*Initialisation des autres champs dans la table des sommets*/
	for(i=0;i<n-1;i++)
	{
		gr->Nodes[i].next_gr = &(gr->Nodes[i+1]);
		gr->Nodes[i].id_W = 1;
		gr->Nodes[i].next_W = &(gr->Nodes[i]);
		gr->Nodes[i].n_sh = &(gr->Nodes[i]);
		gr->Nodes[i].next_sh = &(gr->Nodes[i]);
	}
	gr->Nodes[n-1].next_gr = &(gr->Nodes[0]);
	gr->Nodes[n-1].id_W = 1;
	gr->Nodes[n-1].next_W = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].n_sh = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].next_sh = &(gr->Nodes[n-1]);

	return true;
}



bool copy_Graph(I_Graph *gr,I_Graph *Ogr)
{
    long i,n,m;

    n = Ogr->n_Nodes;
    m = Ogr->m_Edges;

    simpleEdge *e_list;
    e_list = new simpleEdge[m];

    for(i=0;i<m;i++)
    {
        e_list[i].node1 = Ogr->Edges[i].back->adjac->id;
        e_list[i].node2 = Ogr->Edges[i].adjac->id;
        e_list[i].cap = Ogr->Edges[i].cap;
    }

    if(!InitGraph_from_list(e_list,n,m,gr))
    {
        cout << "ERROR:copy_graph(). Erreur pendant la recopie du graphe." << endl;
        return false;
    }

    delete []e_list;

    for(i=0;i<n;i++)
    {
        gr->mark_tab[i] = Ogr->mark_tab[i];
        gr->Nodes[i].terminal = Ogr->Nodes[i].terminal;
    }

    for(int i=0; i<m; i++)
    {
        gr->Edges[i].X = Ogr->Edges[i].X;
        gr->Edges[i].back->X = Ogr->Edges[i].back->X;
        gr->Edges[i].cap = Ogr->Edges[i].cap;
        gr->Edges[i].back->cap = Ogr->Edges[i].back->cap;
    }

//    PrintGraph(Ogr);
//    cout << endl << endl << endl << "AND HIS COPY:" << endl << endl << endl;
//    PrintGraph(gr);

    return true;
}



void PrintGraph(I_Graph *g)
{
	long j,m;

	m = g->m_Edges;

	cout << g->n_Nodes << " " << g->m_Edges << endl;
	for(j=0;j<m;j++)
	{
		cout << g->Edges[m+j].adjac->id << " " << g->Edges[j].adjac->id;
		cout << "    " << g->Edges[j].cap << " " << g->Edges[j].num << endl;
	}
}



/*void PrintGraph_Table(I_Graph *g)
{
	long i,m,n;

	n = g->n_Nodes;
	m = g->m_Edges;

	printf("%ld nodes, %ld edges.\n",g->n_Nodes,g->m_Edges);

	puts("Nodes");
	for(i=0;i<n;i++)
	{
		printf("%ld  %ld  %ld  %ld  ",g->Nodes[i].id,g->Nodes[i].next_gr->id,g->Nodes[i].n_sh->id,g->Nodes[i].next_sh->id);

		if(g->Nodes[i].first_edge != NIL_E)
			printf("E%ld\n",g->Nodes[i].first_edge->num);
		else
			printf("%s\n","-");

	}
}


void PrintReducedGraph(I_Graph *g,FILE *fout)
{
	long first,node_cour,n_contr,v;
	I_Edge *e;
	long i;

	bool *mark_tab = new bool[g->n_Nodes];

	for(i=0;i<g->n_Nodes;i++)
	{
		mark_tab[i] = False;
	}

	first = 1;
	node_cour = 1;
	do
	{
		//Parcours de tous les sommets qui sont contractés sur node_cour
		n_contr = node_cour;

		do
		{
			e = g->Nodes[n_contr-1].first_edge;
			while(e != NIL_E)
			{
				v = e->adjac->id;

				if( (mark_tab[v-1] == False) && (g->Nodes[v-1].n_sh->id != g->Nodes[n_contr-1].n_sh->id) )
				{
					fprintf(fout,"e == %ld %ld  %f\n",g->Nodes[n_contr-1].n_sh->id,g->Nodes[v-1].n_sh->id,e->X);
				}

				e = e->next;
			}

			mark_tab[n_contr - 1] = True;

			n_contr = g->Nodes[n_contr-1].next_sh->id;
		}
		while(n_contr != node_cour);

		node_cour = g->Nodes[node_cour-1].next_gr->id;
	}
	while(node_cour != first);
}*/


/*Contraction des sommets terminaux avec des non-terminaux*/
/*Pour chaque terminal t, on contracte l'arete ut, ssi    */
/*x(delta(u)\t) = x(ut). Cela inclut le cas ou x(ut) = x(uv)*/
void operationTheta(I_Graph *gr,int k)
{
	I_Node *nptr;
	b_Edge *ecour,*ecour2,*delta_u,*delta_t;
	b_Node *W;
	int u,t;
	bool contraction;

	//On commence par contracter les non terminaux
	//les eventuels sommets pendants.
	nptr = gr->Nodes[0].n_sh;
	do
	{
		//Il faut le faire car on peut avoir par exemple nptr->id = 12
		//qu'on contracte avec 6. Comme on insere les sommet dans l'ordre croissant
		//on va contracter l'ensemble (6,12). A l'iteration suivante, on reste sur
		//nptr. Or 12 a disparu et a ete remplace par 6. Donc il faut faire l'instruction
		//suivante pour eviter le probleme.
		nptr = nptr->n_sh;
		contraction = false;

		if(nptr->terminal == false)
		{
			delta_u = get_delta_v(gr,nptr->id,NO_0_EDGES,CAP_VALUE_X);

			//On verifie si le sommet est pendant
			//au lieu de regarder si |delta(W)| = 1,
			//on regarde plutot si le nombre de voisin
			//est egal a 1

			//cout << "NPTR == " << nptr->id << endl;

			int voisinPrec;

			if(delta_u != NIL_BE)
			{
				ecour = delta_u;
				voisinPrec = ecour->node2;

				while((ecour != NIL_BE) && (ecour->node2 == voisinPrec))
				{
					ecour = ecour->next;
				}

				//cout << "Voisin == " << voisinPrec << endl;

				//Il n'y a qu'un seul voisin.
				//Dans ce cas, on contracte u et le voisin
				if(ecour == NIL_BE)
				{
					W = NIL_BN;
					Insere_Node(&W,nptr->id);
					Insere_Node(&W,voisinPrec);
					contraction = true;
				}
			}

			if(contraction)
			{
				/*Print_b_Node_Set(W,stdout);
				cout << endl;*/

				contract_set_w(gr,W);
				Delete_b_Node_Set(&W);
			}

			Delete_b_Edge_Set(&delta_u);
		}


		if(contraction == false)
			nptr = nptr->next_gr;

	}while((contraction == true) ||
		((contraction == false) && (nptr != gr->Nodes[0].n_sh)));

	//Contraction avec les sommets terminaux

	nptr = gr->Nodes[0].n_sh;
	do
	{
		//Il faut le faire car on peut avoir par exemple nptr->id = 12
		//qu'on contracte avec 6. Comme on insere les sommet dans l'ordre croissant
		//on va contracter l'ensemble (6,12). A l'iteration suivante, on reste sur
		//nptr. Or 12 a disparu et a ete remplace par 6. Donc il faut faire l'instruction
		//suivante pour eviter le probleme.
		nptr = nptr->n_sh;
		contraction = false;

		//On regarde s'il s'agit d'un terminal
		if(nptr->terminal == true)
		{
			t = nptr->id;

			//On insere le noeud en precisant
			//que c'est un terminal
			delta_t = get_delta_v(gr,t,NO_0_EDGES,CAP_VALUE_X);

			ecour = delta_t;
			while(ecour != NIL_BE)
			{
				u = ecour->node2;

				if(gr->Nodes[u-1].terminal == false)
				{
					delta_u = get_delta_v(gr,u,NO_0_EDGES,CAP_VALUE_X);

					//Calcul de x(delta(u))
					double sommeX = 0;
					ecour2 = delta_u;
					while(ecour2 != NIL_BE)
					{
						sommeX += ecour2->cap;

						ecour2 = ecour2->next;
					}

					Delete_b_Edge_Set(&delta_u);

					//Si x(delta(u)\ut) = x(ut)
					//On contracte l'arete
					/*cout << "Operation t = " << nptr->id;
					cout << " " << u << endl;
					cout << "SommeX == " << sommeX << endl;
					cout << "Ecour = " << ecour->cap << endl;*/

					//On calcule x([u,t]): il peut y avoir des
					//aretes paralleles
					b_Edge *delta_u_t = get_edge_u_v(gr,t,u,NO_0_EDGES,CAP_VALUE_X);
					double x_u_t = 0;

					ecour2 = delta_u_t;
					while(ecour2 != NIL_BE)
					{
						x_u_t += ecour2->cap;

						ecour2 = ecour2->next;
					}

					Delete_b_Edge_Set(&delta_u_t);

					/*cout << "U == " << u << " T == " << t <<
					endl;
					cout << "Xut == " << x_u_t << " " <<
					sommeX << endl;*/

					if(2*x_u_t >= (sommeX-EPSILON)/*2*ecour->cap >= (sommeX-EPSILON)*/)
					{
						W = NIL_BN;
						Insere_Node(&W,t);
						Insere_Node(&W,u);
						contraction = true;

						/*Print_b_Node_Set(W,stdout);
						cout << endl;*/
					}
				}

				ecour = ecour->next;
			}

			if(contraction)
			{
				contract_set_w(gr,W);
				Delete_b_Node_Set(&W);
			}

			Delete_b_Edge_Set(&delta_t);

		}

		if(contraction == false)
			nptr = nptr->next_gr;



	}while((contraction == true) ||
		((contraction == false) && (nptr != gr->Nodes[0].n_sh)));
}


/*Permet de sauvegarder la contraction courante*/
/*dans les champs annexes*/
void save_contraction_info(I_Graph *gr,int niveau)
{
	long n,i;

	n = gr->n_Nodes;

	if(niveau == 1)
	{
		for(i=0;i<n;i++)
		{
			gr->Nodes[i].a_next_gr = gr->Nodes[i].next_gr;
			gr->Nodes[i].a_n_sh = gr->Nodes[i].n_sh;
			gr->Nodes[i].a_next_sh = gr->Nodes[i].next_sh;
		}
	}
	else if(niveau == 2)
	{
		for(i=0;i<n;i++)
		{
			gr->Nodes[i].b_next_gr = gr->Nodes[i].next_gr;
			gr->Nodes[i].b_n_sh = gr->Nodes[i].n_sh;
			gr->Nodes[i].b_next_sh = gr->Nodes[i].next_sh;
		}
	}
}


/*Permet de réinitialiser la table des contractions*/
/*en rappelant la dernière contraction sauvegardée*/
void recall_contraction_info(I_Graph *gr,int niveau)
{
	long n,i;

	n = gr->n_Nodes;

	if(niveau == 1)
	{
		for(i=0;i<n;i++)
		{
			gr->Nodes[i].next_gr = gr->Nodes[i].a_next_gr;
			gr->Nodes[i].n_sh = gr->Nodes[i].a_n_sh;
			gr->Nodes[i].next_sh = gr->Nodes[i].a_next_sh;
		}
	}
	else if(niveau == 2)
	{
		for(i=0;i<n;i++)
		{
			gr->Nodes[i].next_gr = gr->Nodes[i].b_next_gr;
			gr->Nodes[i].n_sh = gr->Nodes[i].b_n_sh;
			gr->Nodes[i].next_sh = gr->Nodes[i].b_next_sh;
		}
	}
}


/*Retourne Vrai si l'entier z est pair et Faux sinon*/
Bool is_pair(long z)
{
	ldiv_t parite;

	parite = ldiv(z,2);

	if(parite.rem == 1)
		return false;  /// Instead of defined False
	else
		return true;  /// Instead of defined True
}


long sup_part(long num,long denom)
{
	ldiv_t parite;

	parite = ldiv(num,denom);

	return parite.quot + parite.rem;
}


Bool admissible(I_Graph *gr_dem,b_Node *W)
{
	b_Edge *delta_w_dem;

	delta_w_dem = get_delta_w(gr_dem,W,NO_0_EDGES,CAP_VALUE_1);

	Bool res = check_delta_w_card(delta_w_dem,1);

	Delete_b_Edge_Set(&delta_w_dem);

	return(res);
}


Bool admissible(I_Graph *gr_dem,b_Node *partition[],int nbElem)
{
	int i=0;

	while((i < nbElem) && (admissible(gr_dem,partition[i])))
		i++;

	return(i >= nbElem);
}


/*Gives the list of the edges of delta(W).*/
/*The method is similar to the method used in get_edge_list().*/
b_Edge *get_delta_w(I_Graph *g,b_Node *w,int edge_flag,int cap_flag)
{
	b_Node *b_cour;
	long k,nod1,nod2;
	double cap;
	I_Edge *eptr;
	long n_u,n_v,i;
	long num;

	b_Edge *delta_w,*edge_cour;

	if(w == NIL_BN)
	{
		return NIL_BE;
	}

	/*Marking up all the nodes of the set W*/
	if(g->n_Nodes > N_TAB)
	{
		cout << "N nodes superieur à N_TAB" << endl;
		cout << "N = " << g->n_Nodes << " N_TAB = " << N_TAB << endl;

		exit(-1);
	}


	for(i=0;i<g->n_Nodes;i++)
	{
		g->mark_tab[i] = NODE_NOT_MARKED;
	}


	b_cour = w;
	while(b_cour != NIL_BN)
	{
		i = b_cour->id-1;

		g->mark_tab[i] = NODE_MARKED;

		b_cour = b_cour->next;
	}

	/*Construction of the list of edges*/
	/*if(g->m_Edges > M_TAB)
	{
		puts("M edges supérieur à M_TAB");

		exit(-1);
	}*/



	/*We mark 0 on a node to say that this node is not*/
	/*in W. If we mark a non zero value on the node(1 or 2),*/
	/*it means that the node is in W. To build the list E(W),*/
	/*we follow the incidence list of each node of W.*/
	/*When we finish traversing the incidence list of*/
	/*a node u, we mark 2 in mark_tab to say that we*/
	/*have already traverse the incidence list of u.*/
	/*Then when we are traversing the incidence list*/
	/*of a node v(let vk be an edge of this list),we*/
	/*check if k is marked with 2 or not.If k is marked,*/
	/*it means that we have already traverse the incidence*/
	/*list of k, and we have already taken the edge kv.*/
	/*Then we don't take the edge vk.*/
	/*If the node k is not marked with 2, then we can*/
	/*take the edge vk.*/


	/*puts("xxxx");
	for(i=0;i<2*g->m_Edges;i++)
	{
		printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
		//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);

	}*/


	k = 0;

	delta_w = NIL_BE;

	switch(edge_flag)
	{
		case NO_0_EDGES:	/*Suppression des arÃštes nulles: opÃ©ration theta_0*/
			b_cour = w;

			while(b_cour != NIL_BN)
			{
				n_u = b_cour->id;
				do
				{
					eptr = g->Nodes[n_u-1].first_edge;

					while(eptr != NIL_E)
					{
						n_v = eptr->adjac->id;

						if( 	(g->Nodes[n_v-1].n_sh->id != g->Nodes[n_u-1].n_sh->id) 		/*Si u et v ne sont pas contractés ensembles*/
							&&
							(g->mark_tab[g->Nodes[n_v-1].n_sh->id - 1] == NODE_NOT_MARKED) 	/*Si le sommet représentant de v n'est pas marqué*/
							&&
							(eptr->X != 0.0) )						/*Si l'arète est non nulle*/

						{
							nod1 = g->Nodes[n_u-1].n_sh->id;
							nod2 = g->Nodes[n_v-1].n_sh->id;
							num = eptr->num;

							switch(cap_flag)
							{
								case CAP_VALUE_CAP: 	cap = eptr->cap;
											break;

								case CAP_VALUE_X:	cap = eptr->X;
											break;

								case CAP_VALUE_1:	cap = 1.0;
											break;

								case CAP_VALUE_0:	cap = 0.0;
											break;

								case CAP_VALUE_0_FOR_FRAC_EDGE:
											if((1.0 - eptr->X) > 0)	/*Si arÃšte fractionnaire*/
											{
												cap = 0.0;
											}
											else
											{
												cap = eptr->X;
											}
											break;

								case CAP_VALUE_M_FOR_FRAC_EDGE:
											if((1.0 - eptr->X) > 0.0)	/*Si arÃšte fractionnaire*/
											{
												cap = (double)g->m_Edges;
											}
											else
											{
												cap = eptr->X;
											}
											break;


								default:		/*Par défaut on prendra la valeur X*/
											puts("Erreur dans les FLAGS de GET_DELTA_W");
											exit(-1);
											cap = eptr->X;
											break;
							}

							//edge_cour = (b_Edge *)malloc(sizeof(b_Edge));
							edge_cour = new b_Edge;

							edge_cour->node1 = nod1;
							edge_cour->node2 = nod2;
							edge_cour->cap = cap;
							edge_cour->num = num;

							edge_cour->next = delta_w;
							delta_w = edge_cour;


							/*edge_tab[k].node1 = nod1;
							edge_tab[k].node2 = nod2;
							edge_tab[k].cap = cap;

							k++;*/
						}

						eptr = eptr->next;
					}

					n_u = g->Nodes[n_u-1].next_sh->id;

				}while(n_u != b_cour->id);

				/*When we finished the contraction list of b_cour*/
				/*we must set that nodes of the contraction list*/
				/*have all been traverse. So we shouldn't take their*/
				/*incident edges again.*/
				/*mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;*/



				b_cour = b_cour->next;
			}

			break;


		case ALL_EDGES:
			b_cour = w;

			while(b_cour != NIL_BN)
			{
				n_u = b_cour->id;

				do
				{
					eptr = g->Nodes[n_u-1].first_edge;

					while(eptr != NIL_E)
					{
						n_v = eptr->adjac->id;

						if( 	(g->Nodes[n_v-1].n_sh->id != g->Nodes[n_u-1].n_sh->id) 		/*Si u et v ne sont pas contractés ensembles*/
							&&
							(g->mark_tab[g->Nodes[n_v-1].n_sh->id - 1] == NODE_NOT_MARKED) )	/*Si le sommet représentant de v n'est pas marqué*/
						{
							nod1 = g->Nodes[n_u-1].n_sh->id;
							nod2 = g->Nodes[n_v-1].n_sh->id;
							num = eptr->num;

							switch(cap_flag)
							{
								case CAP_VALUE_CAP: 	cap = eptr->cap;
											break;

								case CAP_VALUE_X:	cap = eptr->X;
											break;

								case CAP_VALUE_1:	cap = 1.0;
											break;

								case CAP_VALUE_0:	cap = 0.0;
											break;

								case CAP_VALUE_0_FOR_FRAC_EDGE:
											if((1.0 - eptr->X) > 0.0)	/*Si arÃšte fractionnaire*/
											{
												cap = 0.0;
											}
											else
											{
												cap = eptr->X;
											}
											break;

								case CAP_VALUE_M_FOR_FRAC_EDGE:
											if((1.0 - eptr->X) > 0.0)	/*Si arÃšte fractionnaire*/
											{
												cap = (double)g->m_Edges;
											}
											else
											{
												cap = eptr->X;
											}
											break;


								default:		/*Par défaut on prendra la valeur X*/
											puts("Erreur dans les FLAGS de GET_DELTA_W");
											exit(-1);
											cap = eptr->X;
											break;
							}

							//edge_cour = (b_Edge *)malloc(sizeof(b_Edge));
							edge_cour = new b_Edge;

							edge_cour->node1 = nod1;
							edge_cour->node2 = nod2;
							edge_cour->cap = cap;
							edge_cour->num = num;

							edge_cour->next = delta_w;
							delta_w = edge_cour;


							/*edge_tab[k].node1 = nod1;
							edge_tab[k].node2 = nod2;
							edge_tab[k].cap = cap;

							k++;*/
						}

						eptr = eptr->next;
					}

					n_u = g->Nodes[n_u-1].next_sh->id;

				}while(n_u != b_cour->id);

				/*When we finished the contraction list of b_cour*/
				/*we must set that nodes of the contraction list*/
				/*have all been traverse. So we shouldn't take their*/
				/*incident edges again.*/
				/*mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;*/



				b_cour = b_cour->next;
			}

			break;

	}

	return delta_w;
}


b_Edge *get_edge_u_v(I_Graph *g,long u,long v,int edge_flag,int cap_flag)
{
	b_Edge *becour,*edge_u_v,*be;
	b_Edge *delta_u;

	if(v > g->n_Nodes)
	{
		fprintf(stdout,"%s","ERROR.Get_delta_v. I_Node v is out of bound.");
		return NULL;
	}

	/*Marking up all the nodes of the set W*/
	if(g->n_Nodes > N_TAB)
	{
		cout << "N nodes superieur à N_TAB" << endl;
		cout << "N = " << g->n_Nodes << " N_TAB = " << N_TAB << endl;


		exit(-1);
	}



	/*Construction of the list of edges*/
	/*if(g->m_Edges > M_TAB)
	{
		puts("M edges supérieur à M_TAB");

		exit(-1);
	}*/


	delta_u = get_delta_v(g,u,edge_flag,cap_flag);

	edge_u_v = NIL_BE;
	becour = delta_u;
	while(becour != NIL_BE)
	{
		if(becour->node2 == v)
		{
			be = new b_Edge;
			be->node1 = u;
			be->node2 = v;
			be->cap = becour->cap;
			be->num = becour->num;
			be->next = NIL_BE;

			Insere_Edge(&edge_u_v,be);
		}

		becour = becour->next;
	}

	Delete_b_Edge_Set(&delta_u);

	return edge_u_v;
}



/*Contract the set of nodes given by f_w*/
Bool contract_set_w(I_Graph *gr,b_Node *f_w)
{
	b_Node *b_cour;
	I_Node *n_cour,*n_first,*n_aux;
	long k,i;

	/*We don't contract a set contains less than 2 nodes*/
	if((f_w == NULL) || ((f_w != NULL) && (f_w->next == NULL)))
		return false;  /// Instead of defined False

	/*We start the contraction by creating*/
	/*the list of the node with the field next_sh*/
	/*and by updating the representant of each node*/
	/*of the set. This representant will be the first*/
	/*node of the set.*/

	n_first = &(gr->Nodes[f_w->id - 1]);
	b_cour = f_w->next;
	while(b_cour != NIL_BN)
	{
		k = b_cour->id - 1;
		n_cour = &(gr->Nodes[k]);

		/*Update of the representant of the nodes*/
		/*that have been already contracted in the node n_cour.*/
		/*There is no matter if the node n_cour hasn't been contracted yet.*/
		n_aux = n_cour;
		do
		{
			gr->Nodes[n_aux->id-1].n_sh = n_first;

			/*n_aux->n_sh = n_first;*/

			n_aux = n_aux->next_sh;
		}while(n_aux->id != n_cour->id);

		/*Merging the next_sh list of n_cour to the list of n_first*/
		cyclic_set_merge(n_first,n_cour);


		/*Update of the shrunk graph list*/
		/*We search the previous node of b_cour in the shunk list*/
		i = 0;	/*We suppose that the node 1 is the first of the shrunk graph*/
		while(gr->Nodes[i].next_gr->id != b_cour->id)
		{
			i = gr->Nodes[i].next_gr->id - 1;
		}

		gr->Nodes[i].next_gr = gr->Nodes[k].next_gr;


		/*************/
		b_cour = b_cour->next;
	}

	return true;  /// Instead of defined True
}


/*On cherche ici des triangles admissible en configuration SP*/
void get_sp_partition_chaine_steiner(I_Graph *gr,b_Node **chaine,long *chaine_sz,long *nb_chaine)
{
	long i;
	long u_cour = 0,u = 0,v = 0;
	b_Node *b_cour;
	b_Edge *delta_v,*be_cour;
	I_Node *nptr,*new_node;
	long n,m,nb_voisin,prec;
	Bool nouveau_sommet = true;  /// Instead of defined True
	short *marquage;
	b_Node *dernier;


	n = gr->n_Nodes;
	m = gr->m_Edges;

	marquage = new short[n];

	for(i=0;i<n;i++)
	{
		gr->mark_tab[i] = 0;
	}

	/*Recherche des motifs*/
	/*On met le résultat dans chaine[i]*/
	/*on cherche des motifs de taille 2*/
	(*nb_chaine) = 0;
	nptr = &(gr->Nodes[0]);
	do
	{
		for(i=0;i<n;i++)
		{
			marquage[i] = 0;
		}

		delta_v = get_delta_v(gr,nptr->id,NO_0_EDGES,CAP_VALUE_X);
		be_cour = delta_v;
		while(be_cour != NIL_BE)
		{
			if(gr->Nodes[be_cour->node2-1].terminal == true)
				marquage[be_cour->node2-1]++;

			be_cour = be_cour->next;
		}

		Delete_b_Edge_Set(&delta_v);


		/*Comptage du nombre de voisin par arête double*/
		nb_voisin = 0;
		for(i=0;i<n;i++)
		{
			if((marquage[i] >= 2) && (gr->mark_tab[i] == 0))
			{
				nb_voisin++;
				u_cour = i+1;
			}
		}

		if(nb_voisin == 1)
		{
			/*U contient le sommet en question*/
			chaine[*nb_chaine] = new b_Node;
			chaine[*nb_chaine]->id = nptr->id;

			chaine[*nb_chaine]->next = new b_Node;
			chaine[*nb_chaine]->next->id = u_cour;

			chaine[*nb_chaine]->next->next = NIL_BN;

			dernier = chaine[*nb_chaine]->next;

			gr->mark_tab[nptr->id-1] = 1;
			gr->mark_tab[u_cour-1] = 1;

			chaine_sz[*nb_chaine] = 2;

			prec = nptr->id;
			do
			{
				delta_v = get_delta_v(gr,u_cour,NO_0_EDGES,CAP_VALUE_X);

				for(i=0;i<n;i++)
				{
					marquage[i] = 0;
				}

				be_cour = delta_v;
				while(be_cour != NIL_BE)
				{
					if(gr->Nodes[be_cour->node2-1].terminal == true)
						marquage[be_cour->node2-1]++;

					be_cour = be_cour->next;
				}

				Delete_b_Edge_Set(&delta_v);

				/*Comptage du nombre de voisin par arête double*/
				nb_voisin = 0;
				for(i=0;i<n;i++)
				{
					if((marquage[i] >= 2) && (gr->mark_tab[i] == 0))
					{
						nb_voisin++;
						if(nb_voisin == 1)
							u = i+1;
						else
							v = i+1;
					}
				}


				if(nb_voisin >= 1)
				{
					/*if(u == prec)
					u = v;*/

					b_cour = new b_Node;
					b_cour->id = u;

					/*Il faut insérer en fin de liste*/
					dernier->next = b_cour;
					b_cour->next = NIL_BN;

					dernier = b_cour;

					/*b_cour->next = chaine[*nb_chaine];
					chaine[*nb_chaine] = b_cour;*/

					gr->mark_tab[u-1] = 1;

					/*Longueur de la chaine: pour calculer les RHS*/
					chaine_sz[*nb_chaine]++;
				}

				prec = u_cour;
				u_cour = u;
			}
			while(nb_voisin != 0);


			/*puts("Chaine == ");
			Print_b_Node_Set(chaine[*nb_chaine],stdout);*/
			(*nb_chaine)++;
		}

		/*Recherche d'un nouveau sommet*/
		new_node = nptr;
		do
		{
			new_node = new_node->next_gr;
		}
		while((new_node != gr->Nodes[0].n_sh) && ( (gr->mark_tab[new_node->id-1] != 0) ||
		(new_node->terminal == false)) );


		/*nptr = nptr->next_gr;*/
		if(new_node != gr->Nodes[0].n_sh)
		{
			nouveau_sommet = true;  /// Instead of defined True

			nptr = new_node;
		}
		else
		{
			nouveau_sommet = false;  /// Instead of defined False
		}


	}while(nouveau_sommet == true);   /// Instead of defined True

	delete[] marquage;

	/*puts("CH == ");
	for(i=0;i<(*nb_chaine);i++)
	{
		Print_b_Node_Set(chaine[i],stdout);
	}*/
}


Bool separation_sp_partition_2(I_Graph *gr,int k_ordre,b_Node *chaine[],long *chaine_sz,long
nb_chaine,long ***sp_part_list,long *sp_p,long *nb_sp_part,long *rhs)
{
	b_Edge *sp_part;
	short *mark_edge;
	long i,j,k,n,m;
	b_Edge **delta_Vi;
	b_Edge *be_cour;

	b_Node *b_cour,*b_cour2;
	b_Edge *nouveau;
	I_Node *nptr;
	long *renum1,*renum2;
	Bool trouver;
	long **pcc;
	double X_D;
	long rhs_sp;
	long n1,n2;
	Bool resultat,is_serie_paralelle;
	b_Edge *sp_mop;
	long **partition;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	(*nb_sp_part) = 0;

	mark_edge = new short[m];

	for(j=0;j<nb_chaine;j++)
	{
		delta_Vi = new b_Edge *[chaine_sz[j]];

		for(i=0;i<n;i++)
		{
			gr->mark_tab[i] = 0;
		}


		/*On commence par marquer les sommets de la chaine*/
		b_cour = chaine[j];
		while(b_cour != NIL_BN)
		{
			gr->mark_tab[b_cour->id-1] = 1;

			b_cour =b_cour->next;
		}

		/*Calcul des delta_Vi*/
		k = 0;
		b_cour = chaine[j];
		while(b_cour != NIL_BN)
		{
			delta_Vi[k] = get_delta_v(gr,b_cour->id,ALL_EDGES,CAP_VALUE_X);

			k++;
			b_cour = b_cour->next;
		}


		/*Il faut vérifier que le graphe support est série-paralelle*/
		if(chaine_sz[j] >= 3)
		{

			is_serie_paralelle = true;  /// Instead of defined True
			b_cour = chaine[j];
			k = 0;
			while((k < chaine_sz[j]-1) && (is_serie_paralelle))
			{
				b_cour2 = b_cour->next->next;
				while((b_cour2 != NIL_BN) && (is_serie_paralelle))
				{
					be_cour = delta_Vi[k];
					while((be_cour != NIL_BE) && ((be_cour->cap == 0) || ( (be_cour->cap > 0) && (be_cour->node2
					!= b_cour2->id)) ) )
					{
						be_cour = be_cour->next;
					}

					if(be_cour != NIL_BE)
					{
						is_serie_paralelle = false;  /// Instead of defined False
					}

					b_cour2 = b_cour2->next;
				}

				b_cour = b_cour->next;
				k++;
			}

		}
		else
		{
			is_serie_paralelle = true;  /// Instead of defined True
		}


		if(is_serie_paralelle == true)  /// Instead of defined True
		{
			/*Calcul de la liste des arêtes*/
			for(k=0;k<m;k++)
			{
				mark_edge[k] = 0;
			}

			sp_part = NIL_BE;
			for(k=0;k<chaine_sz[j];k++)
			{
				be_cour = delta_Vi[k];
				while(be_cour != NIL_BE)
				{
					if(mark_edge[be_cour->num-1] == 0)
					{
						//nouveau = (b_Edge *)malloc(sizeof(b_Edge));
						nouveau = new b_Edge;

						nouveau->node1 = be_cour->node1;
						nouveau->node2 = be_cour->node2;
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;

						nouveau->next = sp_part;
						sp_part = nouveau;

						mark_edge[be_cour->num-1] = 1;
					}

					be_cour = be_cour->next;
				}
			}

			/*suppression de la liste des arêtes*/
			/*for(i=0;i<chaine_sz[j];i++)
			{
				Delete_b_Edge_Set(&delta_Vi[i]);
			}
			free(delta_Vi);*/

			/*Calcul des coef*/
			/*p >= 4. Donc il y aura des coef 2*/
			/*p == 3. Tous les coef sont à 1*/
			if(chaine_sz[j] >= 3)
			{
				/*Renumérotation de la liste d'arête*/
				renum1 = new long[n];
				renum2 = new long[chaine_sz[j]+1];
				//renum1 = (long *)malloc(n*sizeof(long));
				//renum2 = (long *)malloc((chaine_sz[j]+1)*sizeof(long));


				b_cour = chaine[j];
				k = 0;
				while(b_cour != NIL_BN)
				{
					k++;
					renum1[b_cour->id-1] = k;
					renum2[k-1] = b_cour->id;

					b_cour = b_cour->next;
				}


				/*MAJ de la numérotation de tous les sommets de W*/
				k++;

				nptr = &(gr->Nodes[0]);
				trouver = false;  /// Instead of defined False
				do
				{
					if(gr->mark_tab[nptr->id-1] == 0)
					{
						renum1[nptr->id-1] = k;
						renum2[k-1] = nptr->id;
					}

					nptr = nptr->next_gr;

				}while( nptr != &(gr->Nodes[0]) );


				/*Calcul des arêtes faisant partie du MOP: y compris les arêtes à 0*/

				sp_mop = NIL_BE;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					/*On ne prend pas une arête entre les sommets i et i+2 au moins*/
					if( ( (renum1[be_cour->node1-1] != k) && (renum1[be_cour->node2-1] != k) &&
					(abs(renum1[be_cour->node2-1] - renum1[be_cour->node1-1]) <= 1) ) ||
					(renum1[be_cour->node1-1] == k) || (renum1[be_cour->node2-1] == k) )
					{
						//nouveau = (b_Edge *)malloc(sizeof(b_Edge));
						nouveau = new b_Edge;

						nouveau->node1 = renum1[be_cour->node1-1];
						nouveau->node2 = renum1[be_cour->node2-1];
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;

						nouveau->next = sp_mop;
						sp_mop = nouveau;
					}

					be_cour->node1 = renum1[be_cour->node1-1];
					be_cour->node2 = renum1[be_cour->node2-1];

					be_cour = be_cour->next;
				}

				/*puts("MOP == ");
				Print_b_Edge_Set(sp_mop,stdout);*/

				pcc = Floyd(sp_mop,k);


				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					n1 = be_cour->node1;
					n2 = be_cour->node2;

					/*n1 = renum2[be_cour->node1-1];
					n2 = renum2[be_cour->node2-1];*/

					X_D = X_D + pcc[n1-1][n2-1]*be_cour->cap;

					/*On stocke le coef de l'arête dans son champ cap*/
					be_cour->cap = pcc[n1-1][n2-1];

					be_cour = be_cour->next;
				}

				delete[] renum1;
				delete[] renum2;
				//free(renum1);
				//free(renum2);

				Delete_b_Edge_Set(&sp_mop);

				for(i=0;i<k;i++)
				{
					delete[] pcc[i];
					//free(pcc[i]);
				}

				delete[] pcc;
				//free(pcc);
			}
			else
			{
				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					X_D = X_D + be_cour->cap;

					be_cour = be_cour->next;
				}
			}


			/*Calcul du RHS*/
			rhs_sp = sup_part(k_ordre,2)*(chaine_sz[j]+1) - 1;

			//printf("RHS == %li  X_D %f\n",rhs_sp,X_D);

			//Print_b_Edge_Set(sp_part,stdout);

			if(X_D < (rhs_sp - EPSILON))
			{
				(*sp_part_list)[*nb_sp_part] = new long[n];

				for(i=0;i<n;i++)
				{
					(*sp_part_list)[*nb_sp_part][i] = chaine_sz[j]+1;
				}

				i = 0;
				b_cour = chaine[j];
				while(b_cour != NIL_BN)
				{
					nptr = &(gr->Nodes[b_cour->id-1]);
					do
					{
						(*sp_part_list)[*nb_sp_part][nptr->id-1] = i+1;

						nptr = nptr->next_sh;

					}while(nptr != &(gr->Nodes[b_cour->id-1]));

					i++;
					b_cour = b_cour->next;
				}


				Delete_b_Edge_Set(&sp_part);

				//puts("SP-CHAINE violée");
				rhs[*nb_sp_part] = rhs_sp;
				sp_p[*nb_sp_part] = chaine_sz[j]+1;

				(*nb_sp_part)++;
			}
			else
			{
				Delete_b_Edge_Set(&sp_part);

				//puts("SP-CHAINE non violée");
			}
		}

		/*suppression de la liste des arêtes*/
		for(i=0;i<chaine_sz[j];i++)
		{
			Delete_b_Edge_Set(&delta_Vi[i]);
		}
		delete[] delta_Vi;
	}

	if((*nb_sp_part) >= 1)
	{
		resultat = true;  /// Instead of defined True
	}
	else
	{
		resultat = false;  /// Instead of defined False
	}

	delete[] mark_edge;

	return resultat;
}


Bool separation_sp_partition_steiner(I_Graph *gr,int k_ordre,b_Node *chaine[],long *chaine_sz,
long nb_chaine,long ***sp_part_list,long *sp_p,long *nb_sp_part,long *rhs)
{
	b_Edge *sp_part;
	short *mark_edge;
	long i,j,k,n,m;
	b_Edge **delta_Vi;
	b_Edge *be_cour;

	b_Node *b_cour,*b_cour2;
	b_Edge *nouveau;
	I_Node *nptr;
	long *renum1,*renum2;
	Bool trouver;
	long **pcc;
	double X_D;
	long rhs_sp;
	long n1,n2;
	Bool resultat,is_serie_paralelle;
	b_Edge *sp_mop;
	long **partition;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	(*nb_sp_part) = 0;

	mark_edge = new short[m];

	for(j=0;j<nb_chaine;j++)
	{
		delta_Vi = new b_Edge *[chaine_sz[j]];

		for(i=0;i<n;i++)
		{
			gr->mark_tab[i] = 0;
		}


		/*On commence par marquer les sommets de la chaine*/
		b_cour = chaine[j];
		while(b_cour != NIL_BN)
		{
			gr->mark_tab[b_cour->id-1] = 1;

			b_cour =b_cour->next;
		}

		/*Calcul des delta_Vi*/
		k = 0;
		b_cour = chaine[j];
		while(b_cour != NIL_BN)
		{
			delta_Vi[k] = get_delta_v(gr,b_cour->id,ALL_EDGES,CAP_VALUE_X);

			k++;
			b_cour = b_cour->next;
		}


		/*Il faut vérifier que le graphe support est série-paralelle*/
		if(chaine_sz[j] >= 3)
		{

			is_serie_paralelle = true;  /// Instead of defined True
			b_cour = chaine[j];
			k = 0;
			while((k < chaine_sz[j]-1) && (is_serie_paralelle))
			{
				b_cour2 = b_cour->next->next;
				while((b_cour2 != NIL_BN) && (is_serie_paralelle))
				{
					be_cour = delta_Vi[k];
					while((be_cour != NIL_BE) && ((be_cour->cap == 0) || ( (be_cour->cap > 0) && (be_cour->node2
					!= b_cour2->id)) ) )
					{
						be_cour = be_cour->next;
					}

					if(be_cour != NIL_BE)
					{
						is_serie_paralelle = false;   /// Instead of defined True or False
					}

					b_cour2 = b_cour2->next;
				}

				b_cour = b_cour->next;
				k++;
			}

		}
		else
		{
			is_serie_paralelle = true;   /// Instead of defined True or False
		}


		if(is_serie_paralelle == true)   /// Instead of defined True or False
		{
			/*Calcul de la liste des arêtes*/
			for(k=0;k<m;k++)
			{
				mark_edge[k] = 0;
			}

			sp_part = NIL_BE;
			for(k=0;k<chaine_sz[j];k++)
			{
				be_cour = delta_Vi[k];
				while(be_cour != NIL_BE)
				{
					if(mark_edge[be_cour->num-1] == 0)
					{
						//nouveau = (b_Edge *)malloc(sizeof(b_Edge));
						nouveau = new b_Edge;

						nouveau->node1 = be_cour->node1;
						nouveau->node2 = be_cour->node2;
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;

						nouveau->next = sp_part;
						sp_part = nouveau;

						mark_edge[be_cour->num-1] = 1;
					}

					be_cour = be_cour->next;
				}
			}

			/*suppression de la liste des arêtes*/
			/*for(i=0;i<chaine_sz[j];i++)
			{
				Delete_b_Edge_Set(&delta_Vi[i]);
			}
			free(delta_Vi);*/


			/*Calcul des coef*/
			/*p >= 4. Donc il y aura des coef 2*/
			/*p == 3. Tous les coef sont à 1*/
			if(chaine_sz[j] >= 3)
			{
				/*Renumérotation de la liste d'arête*/
				renum1 = new long[n];
				renum2 = new long[chaine_sz[j]+1];
				//renum1 = (long *)malloc(n*sizeof(long));
				//renum2 = (long *)malloc((chaine_sz[j]+1)*sizeof(long));


				b_cour = chaine[j];
				k = 0;
				while(b_cour != NIL_BN)
				{
					k++;
					renum1[b_cour->id-1] = k;
					renum2[k-1] = b_cour->id;

					b_cour = b_cour->next;
				}


				/*MAJ de la numérotation de tous les sommets de W*/
				k++;

				nptr = &(gr->Nodes[0]);
				trouver = false;  /// Instead of defined True or False
				do
				{
					if(gr->mark_tab[nptr->id-1] == 0)
					{
						renum1[nptr->id-1] = k;
						renum2[k-1] = nptr->id;
					}

					nptr = nptr->next_gr;

				}while( nptr != &(gr->Nodes[0]) );


				/*Calcul des arêtes faisant partie du MOP: y compris les arêtes à 0*/

				sp_mop = NIL_BE;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					/*On ne prend pas une arête entre les sommets i et i+2 au moins*/
					if( ( (renum1[be_cour->node1-1] != k) && (renum1[be_cour->node2-1] != k) &&
					(abs(renum1[be_cour->node2-1] - renum1[be_cour->node1-1]) <= 1) ) ||
					(renum1[be_cour->node1-1] == k) || (renum1[be_cour->node2-1] == k) )
					{
						//nouveau = (b_Edge *)malloc(sizeof(b_Edge));
						nouveau = new b_Edge;

						nouveau->node1 = renum1[be_cour->node1-1];
						nouveau->node2 = renum1[be_cour->node2-1];
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;

						nouveau->next = sp_mop;
						sp_mop = nouveau;
					}

					be_cour->node1 = renum1[be_cour->node1-1];
					be_cour->node2 = renum1[be_cour->node2-1];

					be_cour = be_cour->next;
				}

				/*puts("MOP == ");
				Print_b_Edge_Set(sp_mop,stdout);*/

				pcc = Floyd(sp_mop,k);


				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					n1 = be_cour->node1;
					n2 = be_cour->node2;

					/*n1 = renum2[be_cour->node1-1];
					n2 = renum2[be_cour->node2-1];*/

					X_D = X_D + pcc[n1-1][n2-1]*be_cour->cap;

					/*On stocke le coef de l'arête dans son champ cap*/
					be_cour->cap = pcc[n1-1][n2-1];

					be_cour = be_cour->next;
				}

				delete[] renum1;
				delete[] renum2;
				//free(renum1);
				//free(renum2);

				Delete_b_Edge_Set(&sp_mop);

				for(i=0;i<k;i++)
				{
					delete[] pcc[i];
					//free(pcc[i]);
				}
				delete[] pcc;
				//free(pcc);
			}
			else
			{
				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while(be_cour != NIL_BE)
				{
					X_D = X_D + be_cour->cap;

					be_cour = be_cour->next;
				}
			}


			/*Calcul du RHS*/
			rhs_sp = sup_part(k_ordre,2)*(chaine_sz[j]+1) - 1;

			//printf("RHS == %li  X_D %f\n",rhs_sp,X_D);

			//Print_b_Edge_Set(sp_part,stdout);

			if(X_D < (rhs_sp - EPSILON))
			{
				(*sp_part_list)[*nb_sp_part] = new long[n];

				for(i=0;i<n;i++)
				{
					(*sp_part_list)[*nb_sp_part][i] = chaine_sz[j]+1;
				}

				i = 0;
				b_cour = chaine[j];
				while(b_cour != NIL_BN)
				{
					nptr = &(gr->Nodes[b_cour->id-1]);
					do
					{
						(*sp_part_list)[*nb_sp_part][nptr->id-1] = i+1;

						nptr = nptr->next_sh;

					}while(nptr != &(gr->Nodes[b_cour->id-1]));

					i++;
					b_cour = b_cour->next;
				}


				Delete_b_Edge_Set(&sp_part);

				//cout << "SP-CHAINE violée" << endl;
				rhs[*nb_sp_part] = rhs_sp;
				sp_p[*nb_sp_part] = chaine_sz[j]+1;

				(*nb_sp_part)++;
			}
			else
			{
				Delete_b_Edge_Set(&sp_part);

				//puts("SP-CHAINE non violée");
			}
		}

		/*suppression de la liste des arêtes*/
		for(i=0;i<chaine_sz[j];i++)
		{
			Delete_b_Edge_Set(&delta_Vi[i]);
		}

		delete[] delta_Vi;
	}

	if((*nb_sp_part) >= 1)
	{
		resultat = true;
	}
	else
	{
		resultat = false;   /// Instead of defined True or False
	}

	delete[] mark_edge;

	return resultat;
}

/*Séparation des Partition*/
/*Retourne Vrai si Partition trouvée et Faux sinon*/
Bool separation_partition_steiner(I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,long **partition,long *rhs,FILE *sortie_frac_gk)
{
	b_Node *V0;
	Bool *mark_tab2;
	long rhs_p,p_card;
	b_Edge **delta_Vi;
	b_Edge *be_cour;
	float X_D;
	b_Node *b_cour;
	long i,s;
	b_Node *b_vi;
	long n,m;
	Bool resultat = false;   /// Instead of defined True or False
	I_Node *nptr;

	/*puts("Part");*/

	n = gr->n_Nodes;
	m = gr->m_Edges;


	V0 = NIL_BN;

	/*Calacul de V0:toujours dans le graphe réduit*/
	//mark_tab2 = (Bool *)calloc(n,sizeof(Bool));
	mark_tab2 = new Bool[n];

	for(i=0;i<n;i++)
	{
		mark_tab2[i] = false;  /// Instead of defined True or False
	}


	/*Marquage des sommets du cycle*/
	b_cour = frac_cycle;
	while(b_cour != NIL_BN)
	{
		mark_tab2[b_cour->id-1] = true;

		b_cour = b_cour->next;
	}

	/*Parcours du graphe réduit*/
	V0 = NIL_BN;
	nptr = gr->Nodes[0].n_sh;
	do
	{
		if(mark_tab2[nptr->id-1] == false)   /// Instead of defined True or False
		{
			Insere_Node(&V0,nptr->id);
		}

		nptr = nptr->next_gr;

	}while( nptr != gr->Nodes[0].n_sh);

	delete[] mark_tab2;


	/*Calcul des delta(Vi):on prend en compte les arètes nulles*/
	if(V0 != NIL_BN)
	{
		p_card = cycle_sz + 1;

		delta_Vi = new b_Edge *[p_card];
		delta_Vi[0] = get_delta_w(gr,V0,ALL_EDGES,CAP_VALUE_X);
		i = 1;
	}
	else
	{
		p_card = cycle_sz;

		delta_Vi = new b_Edge *[p_card];
		i = 0;
	}

	/*Calcul des delta(Vi):on prend en compte les arètes nulles*/
	b_cour = frac_cycle;
	while(b_cour != NIL_BN)
	{
		delta_Vi[i] = get_delta_v(gr,b_cour->id,ALL_EDGES,CAP_VALUE_X);

		i++;
		b_cour = b_cour->next;
	}


	/*Calcul de X_D*/
	X_D = 0;
	for(s=0;s<p_card;s++)
	{
		be_cour = delta_Vi[s];
		while(be_cour != NIL_BE)
		{
			X_D = X_D + be_cour->cap;

			be_cour = be_cour->next;
		}
	}

	/*Calcul du RHS*/
	int p_card_effectif = cycle_sz;

	rhs_p = k*p_card_effectif;
	rhs_p = sup_part(rhs_p,2);


	/*Vérification de la parité de p:il faut que p soit impair*/
	/*puts("V0 == ");
	Print_b_Node_Set(V0,stdout);
	fprintf(sortie_frac_gk,"P == %ld RHS == %ld  X_D == %f\n",p_card,rhs_p,X_D);*/

	resultat = false;   /// Instead of defined True or False

	if(is_pair(k*p_card_effectif) == false)   /// Instead of defined True or False
	{
		/*Vérification de l'inégalité*/
		if( (X_D/2.0) < rhs_p - EPSILON)
		{
			resultat = true;
			//cout << "Partition violée" << endl;

			(*partition) = new long[n];
			for(i=0;i<n;i++)
			{
				(*partition)[i] = 0;
			}

			i = 0;
			b_cour = frac_cycle;
			while(b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id-1]);
				do
				{
					(*partition)[nptr->id-1] = i+1;

					nptr = nptr->next_sh;

				}while(nptr != &(gr->Nodes[b_cour->id-1]));

				b_cour = b_cour->next;
				i++;
			}


			b_cour = V0;
			while(b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id-1]);
				do
				{
					(*partition)[nptr->id-1] = i+1;

					nptr = nptr->next_sh;

				}while(nptr != &(gr->Nodes[b_cour->id-1]));

				b_cour = b_cour->next;
			}
		}
		else
		{
			//fputs("Partition non-violée\n",sortie_frac_gk);
			resultat = false;   /// Instead of defined True or False
		}
	}


	/*Désallocation de la mémoire*/
	Delete_b_Node_Set(&V0);

	for(i=0;i<p_card;i++)
	{
		Delete_b_Edge_Set(&delta_Vi[i]);
	}

	delete[] delta_Vi;


	if(resultat == true)
	{
		(*rhs) = rhs_p;
	}

	return resultat;
}


/*Gives the list of the edges of delta(v).*/
/*The method is similar to the method used in get_delta_w().*/
b_Edge *get_delta_v(I_Graph *g,long v,int edge_flag,int cap_flag)
{
	long k,nod1,nod2;
	double cap;
	I_Edge *eptr;
	long n_u,n_v,num;


	b_Edge *delta_v,*edge_cour;


	if(v > g->n_Nodes)
	{
		fprintf(stdout,"%s","ERROR.Get_delta_v. I_Node v is out of bound.");
		return NULL;
	}

	/*Marking up all the nodes of the set W*/
	if(g->n_Nodes > N_TAB)
	{
		cout << "N nodes superieur à N_TAB" << endl;
		cout << "N = " << g->n_Nodes << " N_TAB = " << N_TAB << endl;

		exit(-1);
	}



	/*Construction of the list of edges*/
	/*if(g->m_Edges > M_TAB)
	{
		puts("M edges supérieur à M_TAB");

		exit(-1);
	}*/

	/*puts("xxxx");
	for(i=0;i<2*g->m_Edges;i++)
	{
		printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
		//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);

	}*/


	k = 0;

	delta_v = NIL_BE;

	switch(edge_flag)
	{
		case NO_0_EDGES:	/*Suppression des arÃštes nulles: opÃ©ration theta_0*/

			n_u = v;
			do
			{
				eptr = g->Nodes[n_u-1].first_edge;

				while(eptr != NIL_E)
				{
					n_v = eptr->adjac->id;

					if( 	(g->Nodes[n_v-1].n_sh->id != g->Nodes[n_u-1].n_sh->id) 		/*Si u et v ne sont pas contractés ensembles*/
						&&
						(eptr->X != 0.0) )						/*Si l'arète est non nulle*/
					{
						nod1 = g->Nodes[n_u-1].n_sh->id;
						nod2 = g->Nodes[n_v-1].n_sh->id;
						num = eptr->num;

						switch(cap_flag)
						{
							case CAP_VALUE_CAP: 	cap = eptr->cap;
										break;

							case CAP_VALUE_X:	cap = eptr->X;
										break;

							case CAP_VALUE_1:	cap = 1L;
										break;

							case CAP_VALUE_0:	cap = 0L;
										break;

							case CAP_VALUE_0_FOR_FRAC_EDGE:
										if((1L - eptr->X) > 0)	/*Si arÃšte fractionnaire*/
										{
											cap = 0L;
										}
										else
										{
											cap = eptr->X;
										}
										break;

							case CAP_VALUE_M_FOR_FRAC_EDGE:
										if((1L - eptr->X) > 0)	/*Si arÃšte fractionnaire*/
										{
											cap = (double)g->m_Edges;
										}
										else
										{
											cap = eptr->X;
										}
										break;


							default:		/*Par défaut on prendra la valeur X*/
										puts("Erreur dans les FLAGS de GET_DELTA_V");
										exit(-1);
										cap = eptr->X;
										break;
						}

						//edge_cour = (b_Edge *)malloc(sizeof(b_Edge));
						edge_cour = new b_Edge;

						edge_cour->node1 = nod1;
						edge_cour->node2 = nod2;
						edge_cour->cap = cap;
						edge_cour->num = num;

						edge_cour->next = delta_v;
						delta_v = edge_cour;


						/*edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;

						k++;*/
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u-1].next_sh->id;

			}while(n_u != v);

			break;


		case ALL_EDGES:
			n_u = v;
			do
			{
				eptr = g->Nodes[n_u-1].first_edge;

				while(eptr != NIL_E)
				{
					n_v = eptr->adjac->id;

					if( (g->Nodes[n_v-1].n_sh->id != g->Nodes[n_u-1].n_sh->id) )		/*Si u et v ne sont pas contractés ensembles*/
					{
						nod1 = g->Nodes[n_u-1].n_sh->id;
						nod2 = g->Nodes[n_v-1].n_sh->id;
						num = eptr->num;

						switch(cap_flag)
						{
							case CAP_VALUE_CAP: 	cap = eptr->cap;
										break;

							case CAP_VALUE_X:	cap = eptr->X;
										break;

							case CAP_VALUE_1:	cap = 1.0;
										break;

							case CAP_VALUE_0:	cap = 0.0;
										break;

							case CAP_VALUE_0_FOR_FRAC_EDGE:
										if((1.0 - eptr->X) > 0)	/*Si arÃšte fractionnaire*/
										{
											cap = 0.0;
										}
										else
										{
											cap = eptr->X;
										}
										break;

							case CAP_VALUE_M_FOR_FRAC_EDGE:
										if((1.0 - eptr->X) > 0)	/*Si arÃšte fractionnaire*/
										{
											cap = (double)g->m_Edges;
										}
										else
										{
											cap = eptr->X;
										}
										break;


							default:		/*Par défaut on prendra la valeur X*/
										puts("Erreur dans les FLAGS de GET_DELTA_V");
										exit(-1);
										cap = eptr->X;
										break;
							}

						//edge_cour = (b_Edge *)malloc(sizeof(b_Edge));
						edge_cour = new b_Edge;

						edge_cour->node1 = nod1;
						edge_cour->node2 = nod2;
						edge_cour->cap = cap;
						edge_cour->num = num;


						edge_cour->next = delta_v;
						delta_v = edge_cour;


						/*edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;

						k++;*/
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u-1].next_sh->id;

			}while(n_u != v);

			break;

	}

	return delta_v;
}

void Dijkstra(I_Graph *gr,int s,double *&L,int *&P,int *&PredArete,int cap_flag)
{
	int n = gr->n_Nodes;

	L = new double[n];
	P = new int[n];
	PredArete = new int[n];

	/*************************/
	//Initialisation de P et L
	b_Edge *ecour,*delta_v;
	int u;

	for(int i=0;i<n;i++)
	{
		L[i] = INFINI;
		P[i] = -1;
		PredArete[i] = -1;
	}

	L[s-1] = 0;

	delta_v = get_delta_v(gr,s,NO_0_EDGES,cap_flag);

	ecour = delta_v;
	while(ecour != NIL_BE)
	{
		u = ecour->node2;

		L[u-1] = ecour->cap;
		P[u-1] = s;
		PredArete[u-1] = ecour->num;

		ecour = ecour->next;
	}

	Delete_b_Edge_Set(&delta_v);

	vector<int> S;

	int nbNoeudActif;
	int *T = new int[n];

	S.push_back(s);

	for(int i=0;i<n;i++)
		T[i] = i+1;

	T[s-1] = -1;

	nbNoeudActif = n-1;

	/**************************/

	while(nbNoeudActif > 0)
	{
		u = trouverMinimum(L,T,n);
		u++;

		T[u-1] = -1;
		nbNoeudActif--;

		S.push_back(u);

		delta_v = get_delta_v(gr,u,NO_0_EDGES,cap_flag);

		ecour = delta_v;
		while(ecour != NIL_BE)
		{
			//cout << "l1[" << ecour->node2 << "] = " <<
			//L[ecour->node2-1] << " " << L[u-1] << " " << ecour->cap << endl;

			//Si le sommet n'est pas dans T
			if((ecour->cap >= 0) && (T[ecour->node2-1] != -1))
			{
				if(L[ecour->node2-1] > L[u-1] + ecour->cap)
				{
					L[ecour->node2-1] = L[u-1] + ecour->cap;
					P[ecour->node2-1] = u;
					PredArete[ecour->node2-1] = ecour->num;
				}
			}

			ecour = ecour->next;
		}

		Delete_b_Edge_Set(&delta_v);
	}

	delete[] T;
}


int trouverMinimum(double *L,int *T,int n)
{
	int i;
	int imin;
	int u,umin;

	imin = 0;
	while(T[imin] == -1)
		imin++;

	for(i=imin+1;i<n;i++)
	{
		if(T[i] != -1)
		{
			u = T[i];
			umin = T[imin];

			if(L[u-1] < L[umin-1])
				imin = i;
		}
	}

	return imin;
}


/*S�paration des Partition*/
/*Retourne Vrai si Partition trouv�e et Faux sinon*/
Bool separation_double_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t,int *&partition,int &F,int &rhs,FILE *sortie_frac_gk)
{
	int n,m;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	//On commence par initialiser la partition
	partition = new int[n];
	for(int i=0;i<5;i++)
	{
		b_Node *bcour;
		bcour = pi[i];
		while(bcour != NIL_BN)
		{
			partition[bcour->id-1] = i;

			bcour = bcour->next;
		}
	}

	//Parcours des aretes du graphe pour verifier
	//si la partition induit une contrainte violee

	//Si k est pair, alors on cherche une arete F
	if(is_pair(k) == true) /// Instead of defined True or False
	{
		b_Edge *delta_v1 = get_delta_w(gr,pi[2],NO_0_EDGES,CAP_VALUE_X);
		b_Edge *ecour,*emax;

		emax = delta_v1;
		ecour = delta_v1->next;

		while(ecour != NIL_BE)
		{
			if(((partition[ecour->node2-1] == 1) || (partition[ecour->node2-1] == 4)) &&
			(ecour->cap > emax->cap))
				emax = ecour;

			ecour = ecour->next;
		}

		F = emax->num;
		Delete_b_Edge_Set(&delta_v1);

	}
	else
		F = -1;

	//On verifie si la contrainte est violee ou pas
	double X_D = 0;
	double coe;

	for(int i=0;i<m;i++)
	{
		if(gr->Edges[i].num == F)
			coe = 0;
		else if(abs(partition[gr->Edges[i].back->adjac->id-1] -
		partition[gr->Edges[i].adjac->id-1]) > 1)
			coe = 1;
		else if(abs(partition[gr->Edges[i].back->adjac->id-1] -
		partition[gr->Edges[i].adjac->id-1]) == 0)
			coe = 0;
		else if(abs(partition[gr->Edges[i].back->adjac->id-1] -
		partition[gr->Edges[i].adjac->id-1]) == 1)
		{
			if(partition[gr->Edges[i].back->adjac->id-1] == 2 && partition[gr->Edges[i].adjac->id-1] == 3)
				coe = 0;
			else if(partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 2)
				coe = 0;
			else if(partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 4)
				coe = 0;
			else if(partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 3)
				coe = 0;
			else
				coe = 1;
		}

		X_D += coe*gr->Edges[i].X;
	}

	Bool resultat;

	rhs = 3*k;
	rhs = sup_part(rhs,2);

	if(X_D < rhs - EPSILON)
	{
		//cout << "Double Cut violee " << endl;
		resultat = true;    /// Instead of defined True or False
	}
	else
	{
		//cout << "Double Cut non violee " << endl;
		resultat = false;    /// Instead of defined True or False

		delete[] partition;
	}

	return resultat;


}


/*S�paration des Partition*/
/*Retourne Vrai si Partition trouv�e et Faux sinon*/
Bool separation_double_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t,int *&partition,vector<int> &F,int &rhs,FILE *sortie_frac_gk)
{
	int n,m;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	//On commence par initialiser la partition
	partition = new int[n];
	for(int i=0;i<5;i++)
	{
		//cout << "V" << i << " = ";

		b_Node *bcour;
		bcour = pi[i];
		while(bcour != NIL_BN)
		{
			//cout << bcour->id << ",";

			partition[bcour->id-1] = i;
			bcour = bcour->next;
		}

		//cout << endl;
	}

	//cout << endl;

	//Parcours des aretes du graphe pour verifier
	//si la partition induit une contrainte violee

	//On cherche k-1 aretes à mettre dans F
	//Notons que quelque soit la parite de k, il faut k-1 aretes

	b_Edge *delta_v1 = get_delta_w(gr,pi[2],NO_0_EDGES,CAP_VALUE_X);
	b_Edge *ecour;

	//On trie les aretes de delta_v1
	vector<int> FauxNum;
	vector<double> FauxCap;
	FauxNum.clear();
	FauxCap.clear();

	ecour = delta_v1;
	while(ecour != NIL_BE)
	{
		if((partition[ecour->node2-1] == 1) || (partition[ecour->node2-1] == 4))
		{

			FauxNum.push_back(ecour->num);
			FauxCap.push_back(ecour->cap);


			/*cout << "ECOUR = " << ecour->num << endl;
			cout << "ECOUR = " << ecour->cap << endl;*/
		}

		ecour = ecour->next;
	}

	Delete_b_Edge_Set(&delta_v1);

	//cout << "FAUX.SIZE = " << FauxNum.size() << endl;

	//Triage des aretes dans Faux
	for(int i=0;i<(int)FauxNum.size()-1;i++)
	{
		int jMax = i;
		for(int j=i+1;j<(int)FauxNum.size();j++)
		{
			if(FauxCap[j] > FauxCap[jMax])
			{
				jMax = j;
			}
		}

		int tmpNum = FauxNum[i];
		FauxNum[i] = FauxNum[jMax];
		FauxNum[jMax] = tmpNum;

		double tmpCap = FauxCap[i];
		FauxCap[i] = FauxCap[jMax];
		FauxCap[jMax] = tmpCap;
	}

	//On transfert ensuite les k-1 premiers dans F
	//Attention, il faut qu'il soit entre V2 et V1
	//ou entre V2 et V4
	Bool resultat = false;    /// Instead of defined True or False

	if(FauxNum.size() >= k-1)
	{
		F.clear();
		for(int i=0;i<k-1;i++)
		{
			F.push_back(FauxNum[i]);
			//cout << "F CAP = " << FauxCap[i] << endl;
		}

		//On verifie si la contrainte est violee ou pas
		double X_D = 0;
		double coe;

		for(int i=0;i<m;i++)
		{
			int j = 0;
			while(j < F.size() && F[j] != gr->Edges[i].num)
				j++;

			if(j < F.size())
				coe = 0;
			else if(abs(partition[gr->Edges[i].back->adjac->id-1] - partition[gr->Edges[i].adjac->id-1]) > 1)
				coe = 1;
			else if(abs(partition[gr->Edges[i].back->adjac->id-1] - partition[gr->Edges[i].adjac->id-1]) == 0)
				coe = 0;
			else if(abs(partition[gr->Edges[i].back->adjac->id-1] - partition[gr->Edges[i].adjac->id-1]) == 1)
			{
				if(partition[gr->Edges[i].back->adjac->id-1] == 2 && partition[gr->Edges[i].adjac->id-1] == 3)
					coe = 0;
				else if(partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 2)
					coe = 0;
				else if(partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 4)
					coe = 0;
				else if(partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 3)
					coe = 0;
				else
					coe = 1;
			}

			X_D += coe*gr->Edges[i].X;
		}

		rhs = 3*k-F.size();
		rhs = sup_part(rhs,2);

		if(X_D < rhs - EPSILON)
		{
			/*cout << "X_D = " << X_D << " RHS = " << rhs << endl;
			cout << "Double Cut violee " << endl;*/
			resultat = true;    /// Instead of defined True or False
		}
		else
		{
			//cout << "Double Cut non violee " << endl;
			resultat = false;    /// Instead of defined True or False

			delete[] partition;
		}
	}
	else
	{
		delete[] partition;
	}

	return resultat;
}

Bool separation_triple_path_cut(I_Graph *gr,b_Node *pi[],int k,int s,int t1,int t2,int *&partition,vector<int> &F,int &rhs,FILE *sortie_frac_gk)
{
	int n,m;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	//On commence par initialiser la partition
	partition = new int[n];
	for(int i=0;i<6;i++)
	{
		//cout << "V" << i << " = ";

		b_Node *bcour;
		bcour = pi[i];
		while(bcour != NIL_BN)
		{
			//cout << bcour->id << ",";

			partition[bcour->id-1] = i;
			bcour = bcour->next;
		}

		//cout << endl;
	}

	//cout << endl;

	//Parcours des aretes du graphe pour verifier
	//si la partition induit une contrainte violee

	//On cherche k-1 aretes à mettre dans F
	//Notons que quelque soit la parite de k, il faut k-1 aretes

	b_Edge *delta_v4 = get_delta_w(gr,pi[4],NO_0_EDGES,CAP_VALUE_X);
	b_Edge *delta_v5 = get_delta_w(gr,pi[5],NO_0_EDGES,CAP_VALUE_X);
	b_Edge *ecour;

	//On trie les aretes de delta_v4
	vector<int> FauxNum;
	vector<double> FauxCap;
	FauxNum.clear();
	FauxCap.clear();

	ecour = delta_v4;
	while(ecour != NIL_BE)
	{
		if((partition[ecour->node2-1] == 2) || (partition[ecour->node2-1] == 3))
		{
			FauxNum.push_back(ecour->num);
			FauxCap.push_back(ecour->cap);

			/*cout << "ECOUR = " << ecour->num << endl;
			cout << "ECOUR = " << ecour->cap << endl;*/
		}

		ecour = ecour->next;
	}

	ecour = delta_v5;
	while(ecour != NIL_BE)
	{
		if(partition[ecour->node2-1] == 3)
		{
			FauxNum.push_back(ecour->num);
			FauxCap.push_back(ecour->cap);

			/*cout << "ECOUR = " << ecour->num << endl;
			cout << "ECOUR = " << ecour->cap << endl;*/
		}

		ecour = ecour->next;
	}

	Delete_b_Edge_Set(&delta_v4);
	Delete_b_Edge_Set(&delta_v5);

	//cout << "FAUX.SIZE = " << FauxNum.size() << endl;

	//Triage des aretes dans Faux
	for(int i=0;i<(int)FauxNum.size()-1;i++)
	{
		int jMax = i;
		for(int j=i+1;j<(int)FauxNum.size();j++)
		{
			if(FauxCap[j] > FauxCap[jMax])
			{
				jMax = j;
			}
		}

		int tmpNum = FauxNum[i];
		FauxNum[i] = FauxNum[jMax];
		FauxNum[jMax] = tmpNum;

		double tmpCap = FauxCap[i];
		FauxCap[i] = FauxCap[jMax];
		FauxCap[jMax] = tmpCap;
	}

	//On transfert ensuite les k-1 premiers dans F
	//Attention, il faut qu'il soit dans [V2 U V3,V4] ou dans [V3,V5]
	Bool resultat = false;    /// Instead of defined True or False

	if(FauxNum.size() >= k-1)
	{
		int z=0;
		F.clear();
		double X_F = 0;
		while((z < FauxNum.size()) && (F.size()-X_F < 1) && (F.size() < k-1))
		{
			F.push_back(FauxNum[z]);
			X_F += FauxCap[z];
			z++;
		}

		if(is_pair(F.size()) == false)   /// Instead of defined True or False
			F.pop_back();

		//cout << "X_F = " << X_F << " |F| = " << F.size() << endl;


		/*for(int i=0;i<k-1;i++)
		{
			F.push_back(FauxNum[i]);
			//cout << "F CAP = " << FauxCap[i] << endl;
		}*/

		//On verifie si la contrainte est violee ou pas
		double X_D = 0;
		double coe;

		for(int i=0;i<m;i++)
		{
			int j = 0;
			while(j < F.size() && F[j] != gr->Edges[i].num)
				j++;

			if(j < F.size())
				coe = 0;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 4)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 3) )
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 2 && partition[gr->Edges[i].adjac->id-1] == 4)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 2) )
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 5)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 5 && partition[gr->Edges[i].adjac->id-1] == 3) )
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 0 && partition[gr->Edges[i].adjac->id-1] == 2)
				|| (partition[gr->Edges[i].back->adjac->id-1] == 2 && partition[gr->Edges[i].adjac->id-1] == 0) )
				coe = 2;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 0 && partition[gr->Edges[i].adjac->id-1] == 3)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 0) )
				coe = 2;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 1 && partition[gr->Edges[i].adjac->id-1] == 3)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 3 && partition[gr->Edges[i].adjac->id-1] == 1) )
				coe = 2;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 0 && partition[gr->Edges[i].adjac->id-1] == 5)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 5 && partition[gr->Edges[i].adjac->id-1] == 0))
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 0 && partition[gr->Edges[i].adjac->id-1] == 4)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 0))
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 1 && partition[gr->Edges[i].adjac->id-1] == 4)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 1))
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 1 && partition[gr->Edges[i].adjac->id-1] == 5)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 5 && partition[gr->Edges[i].adjac->id-1] == 1) )
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 2 && partition[gr->Edges[i].adjac->id-1] == 5)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 5 && partition[gr->Edges[i].adjac->id-1] == 2) )
				coe = 1;
			else if( (partition[gr->Edges[i].back->adjac->id-1] == 4 && partition[gr->Edges[i].adjac->id-1] == 5)
					|| (partition[gr->Edges[i].back->adjac->id-1] == 5 && partition[gr->Edges[i].adjac->id-1] == 4) )
				coe = 1;
			else
				coe = 0;

			X_D += coe*gr->Edges[i].X;
		}

		rhs = 3*k-F.size();
		rhs = sup_part(rhs,2);

		//cout << "X_D = " << X_D << " RHS = " << rhs << endl;

		if(X_D < rhs - EPSILON)
		{
			//cout << "X_D = " << X_D << " RHS = " << rhs << endl;
			//cout << "Triple path cut violee " << endl;
			resultat = true;    /// Instead of defined True or False
		}
		else
		{
			//cout << "Triple path cut non violee " << endl;
			resultat = false;    /// Instead of defined True or False

			delete[] partition;
		}
	}
	else
	{
		delete[] partition;
	}

	return resultat;
}



/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_w_card(b_Node *w,long a)
{
	b_Node *b_cour;
	long card;

	card = 0;
	b_cour = w;
	while((b_cour != NIL_BN) && (card < a))
	{
		card++;
		b_cour = b_cour->next;
	}

	if(card >= a)
		return true;    /// Instead of defined True or False
	else
		return false;    /// Instead of defined True or False
}




/*Checks if the cardinal of the delta(W) is greater or equal to a.*/
/*You can use this function to check if |delta(W)| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_delta_w_card(b_Edge *d_w,long a)
{
	b_Edge *b_cour;
	long card;

	card = 0;
	b_cour = d_w;
	while((b_cour != NIL_BE) && (card < a))
	{
		card++;
		b_cour = b_cour->next;
	}

	if(card >= a)
		return true;    /// Instead of defined True or False
	else
		return false;    /// Instead of defined True or False
}

/*Makes the fusion of two cyclic linked lists*/
/*The lists must have at least 1 node.*/
Bool cyclic_set_merge(I_Node *f1,I_Node *f2)
{
	I_Node *aux;

	if((f1 == NIL_N) || (f2 == NIL_N))
		return false;    /// Instead of defined True or False

	aux = f2->next_sh;
	f2->next_sh = f1->next_sh;
	f1->next_sh = aux;

	return true;    /// Instead of defined True or False
}

/*Calcule le plus court chemin entre*/
/*tous les couples de sommets dans le graphe G*/
/*On utilise l'algorithme de Floyd*/
long **Floyd(b_Edge *sh_eList,long n)
{
	long **pcc;
	long i,j,k;
	long n1,n2;
	b_Edge *be_cour;

	//printf("N == %d\n",n);


	/*Initialisation du tableau des pcc*/
	//pcc = (long **)malloc(n*sizeof(long *));
	pcc = new long* [n];

	for(i=0;i<n;i++)
	{
		//pcc[i] = (long *)malloc(n*sizeof(long));
		pcc[i] = new long[n];
	}

	/*Initialisation du tableau des pcc*/
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			pcc[i][j] = INFINI;
		}
	}



	/*for(i=0;i<m;i++)
	{
		n1 = sh_eList[i].node1;
		n2 = sh_eList[i].node2;

		pcc[n1-1][n2-1] = 1;
		pcc[n2-1][n1-1] = 1;
	}*/


	be_cour = sh_eList;
	while(be_cour != NIL_BE)
	{
		n1 = be_cour->node1;
		n2 = be_cour->node2;

		pcc[n1-1][n2-1] = 1;
		pcc[n2-1][n1-1] = 1;

		be_cour = be_cour->next;
	}


	for(i=0;i<n;i++)
	{
		pcc[i][i] = 0;
	}

	/*nptr1 = &(gr->Nodes[0]);
	do
	{
		k = nptr1->id-1;

		nptr2 = &(gr->Nodes[0]);
		do
		{
			i = nptr2->id-1;

			nptr3 = &(gr->Nodes[0]);
			do
			{
				j = nptr3->id-1;

				pcc[i][j] = min(pcc[i][j],pcc[i][k] + pcc[k][j]);

				nptr3 = nptr3->next_gr;

			}while(nptr3 != &(gr->Nodes[0]));

			nptr2 = nptr2->next_gr;

		}while(nptr2 != &(gr->Nodes[0]));

		nptr1 = nptr1->next_gr;

	}while(nptr1 != &(gr->Nodes[0]));*/


	for(k=0;k<n;k++)
	{
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
				pcc[i][j] = (long)min(pcc[i][j],pcc[i][k] + pcc[k][j]);
			}
		}
	}


	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			//printf("pcc[%d][%d] == %d\n",i,j,pcc[i][j]);
		}
	}

	/*nptr1 = &(gr->Nodes[0]);
	do
	{
		i = nptr1->id-1;

		nptr2 = &(gr->Nodes[0]);
		do
		{
			j = nptr2->id-1;

			printf("i = %d  j = %d  %d\n",i+1,j+1,pcc[i][j]);

			nptr2 = nptr2->next_gr;

		}while(nptr2 != &(gr->Nodes[0]));

		nptr1 = nptr1->next_gr;

	}while(nptr1 != &(gr->Nodes[0]));*/


	return pcc;
}

simpleEdge *GetReducedGraph_eList(I_Graph *g,long *m_edge)
{
	long first,node_cour,n_contr,v;
	I_Edge *e;
	long i,k;
	simpleEdge *eList,*e_aux;

	eList = new simpleEdge[g->m_Edges];

	for(i=0;i<g->n_Nodes;i++)
	{
		g->mark_tab[i] = false;    /// Instead of defined True or False
	}

	k = 0;
	first = 1;
	node_cour = 1;

	do
	{
		/*Parcours de tous les sommets qui sont contractés sur node_cour*/
		n_contr = node_cour;

		do
		{
			e = g->Nodes[n_contr-1].first_edge;
			while(e != NIL_E)
			{
				v = e->adjac->id;

				/*On ne prend pas les arètes à 0*/
				if( (g->mark_tab[v-1] == false) && (g->Nodes[v-1].n_sh->id != g->Nodes[n_contr-1].n_sh->id)    /// Instead of defined True or False
					&& (e->X != 0.0) )
				{

					eList[k].node1 = g->Nodes[n_contr-1].n_sh->id;
					eList[k].node2 = g->Nodes[v-1].n_sh->id;
					eList[k].cap = e->X;

					k++;
				}

				e = e->next;
			}

			g->mark_tab[n_contr - 1] = true;    /// Instead of defined True or False

			n_contr = g->Nodes[n_contr-1].next_sh->id;
		}
		while(n_contr != node_cour);

		node_cour = g->Nodes[node_cour-1].next_gr->id;
	}
	while(node_cour != first);

	/*Optimisation de l'espace*/
	e_aux = eList;

	//eList = (simpleEdge *)malloc(k*sizeof(simpleEdge));
	eList = new simpleEdge[k];

	for(i=0;i<k;i++)
	{
		eList[i].node1 = e_aux[i].node1;
		eList[i].node2 = e_aux[i].node2;
		eList[i].cap = e_aux[i].cap;
	}

	(*m_edge) = k;

	//free(e_aux);
	delete[] e_aux;

	return eList;
}


bool InitGraph_from_vect(vector<demande> *vdem,long n,long m,I_Graph *gr)
{
	int i,j;
	I_Edge *eptr1,*eptr2;
	long node1,node2;
	double cap;
	I_Node *nptr1,*nptr2;

	if(!AllocateGraph(n,m,gr))
	{
		return false;   /// Instead of defined True or False
	}

	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0]);
	eptr2 = &(gr->Edges[m]);

	for(j=0;j<m;j++)
	{
		node1 = (*vdem)[j].source;
		node2 = (*vdem)[j].dest;
		cap = (*vdem)[j].val;

		--node1;
		--node2;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j+1;
		eptr2->num = j+1;

		if(nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NULL;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if(nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NULL;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;

	}

	for(i=0;i<n;i++)
		gr->Nodes[i].terminal = false;

	/*Initialisation des autres champs dans la table des sommets*/
	for(i=0;i<n-1;i++)
	{
		gr->Nodes[i].next_gr = &(gr->Nodes[i+1]);
		gr->Nodes[i].id_W = 1;
		gr->Nodes[i].next_W = &(gr->Nodes[i]);
		gr->Nodes[i].n_sh = &(gr->Nodes[i]);
		gr->Nodes[i].next_sh = &(gr->Nodes[i]);
	}
	gr->Nodes[n-1].next_gr = &(gr->Nodes[0]);
	gr->Nodes[n-1].id_W = 1;
	gr->Nodes[n-1].next_W = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].n_sh = &(gr->Nodes[n-1]);
	gr->Nodes[n-1].next_sh = &(gr->Nodes[n-1]);

	return true;
}


void DrawGraphViz(I_Graph *g,char *nomfich,int parN,vector<simpleEdge> *label,int nbDem)
{
	int p,q;

	I_Edge *ecour;
	ofstream fich(nomfich);

	fich << "graph G{" << endl;

	/*Fixation de la position des sommets et du label*/
	fich << "node [shape=circle];" << endl;

	if(label != NULL)
	{
		for(int i=0;i<g->n_Nodes;i++)
		{
			char s[3];
			char t[3];

            //On calcule d'abord le numero de la demande à laquelle appartient cet arc
            int r,j;
            int parDG = 0;

            bool trouver = false;
            r = 0;
            while(r<nbDem && !trouver)
            {
                j = 0;
                while(j < label[r].size() && parDG+j+1 != i+1)
                    j++;

                if(j < label[r].size())
                    trouver = true;
                else
                {
                    parDG += label[r].size();
                    r++;
                }
            }

            //A la fin, r et j nous donne la position du sommet dans la table des labels

			/*p = (*label)[i].node1;
			q = (*label)[i].node2; */

            p = label[r][j].node1;
			q = label[r][j].node2;

			if(p <= parN)
			{
				p = p;
				s[0] = '\0';
			}
			else if(p <= 2*parN)
			{
				p = p-parN;
				s[0] = '\'';
				s[1] = '\0';
			}
			else if(p <= 3*parN)
			{
				p = p-2*parN;

				s[0] = '\'';
				s[1] = '\'';
				s[2] = '\0';
			}
			else
			{	p = p-3*parN;
				s[0] = '\0';
			}


			if(q <= parN)
			{
				q = q;
				t[0] = '\0';
			}
			else if(q <= 2*parN)
			{
				q = q-parN;
				t[0] = '\'';
				t[1] = '\0';
			}
			else if(q <= 3*parN)
			{
				q = q-2*parN;

				t[0] = '\'';
				t[1] = '\'';
				t[2] = '\0';
			}
			else
			{
				q = q-3*parN;
				t[0] = '\0';
			}

			fich << g->Nodes[i].id << " [label = \"" <<
			p << s << "," << q << t << "," << r+1 << "\"];" << endl;
		}
	}
	else
	{
		for(int i=0;i<g->n_Nodes;i++)
			fich << g->Nodes[i].id << ";" << endl;
	}

	for(int i=0;i<g->m_Edges;i++)
	{
		if(g->Edges[i].X >= 1)
		{
			fich << g->Edges[i].back->adjac->id << "--" << g->Edges[i].adjac->id << "[id=\"";
			fich << g->Edges[i].num << "\"];" << endl;
		}
		else if(g->Edges[i].X > 0)
		{
			fich << g->Edges[i].back->adjac->id << "--" << g->Edges[i].adjac->id;
			fich << "[id=\"" << g->Edges[i].num << "\"]";
			fich << "[style=\"dashed\"][color=\"red\"];" << endl;
		}
	}

	fich << "}" << endl;

	fich.close();
}

