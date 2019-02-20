#include "I_common.h"
#include <math.h>

#include <ilcplex/ilocplex.h>

ostream &operator<<(ostream &os,const simpleEdge &se)
{
	os << se.node1 << " " << se.node2 << " " << se.cap;
}

/*Print the table containing the tree given in parameter*/
void PrintGHCutTree_Table(Tree *ghct)
{
	long i;

	for(i=0;i<ghct->n_nodes;i++)
	{
		cout << ghct->Tab[i].id << " ";

		if(ghct->Tab[i].lv != NULL)
			cout << ghct->Tab[i].lv->id << " ";
		else cout << "-  ";

		if(ghct->Tab[i].lh != NULL)
			cout << ghct->Tab[i].lh->id << " ";
		else cout << "-  ";

		cout << ghct->Tab[i].cap << endl;
	}
}


void PrintGHCutTree(Tree *ghct)
{
	PrintTree(&(ghct->Tab[0]));
}



void PrintTree(Tree_Node *t)
{
	Tree_Node *t_son;

	if(t != NIL_TN)
	{
		cout << t->id << " ";

		t_son = t->lv;
		while(t_son != NIL_TN)
		{
			PrintTree(t_son);

			t_son = t_son->lh;
		}
	}
}

/*Delete the memory allocated for a b_Node set*/
void Delete_b_Node_Set(b_Node **f_w)
{
	b_Node *b_cour;

	while((*f_w) != NIL_BN)
	{
		b_cour = (*f_w);
		(*f_w) = (*f_w)->next;

		//printf("BCOUR == %d\n",b_cour->id);

		delete b_cour;

		//printf("BCOUR 222 == %d\n",b_cour->id);

		//b_cour = NIL_BN;
	}
}


/*Delete the memory allocated for a b_Node set*/
void Delete_b_Edge_Set(b_Edge **f_w)
{
	b_Edge *b_cour;

	while(*f_w != NIL_BE)
	{
		b_cour = (*f_w);
		(*f_w) = (*f_w)->next;

		delete b_cour;

		b_cour = NIL_BE;
	}
}

/*Delete the memory allocated to the tree*/
void Delete_Tree(Tree **tr)
{
	if((*tr) != NULL)
	{
		delete((*tr)->Tab);
		(*tr)->Tab = NULL;

		Delete_b_Edge_Set(&(*tr)->b_List);

		delete (*tr);
	}
}


/*delete the memory allocated to a global mincut*/
void Delete_Mincut(Mincut *gmcu)
{
	Delete_b_Node_Set(&gmcu->f_w);
	Delete_b_Node_Set(&gmcu->f_w_b);
	delete gmcu;
}

/*Supprime une LCC de la mémoire*/
void Delete_LCC_Node(b_Node **f_w)
{
	b_Node *b_cour;

	while((*f_w)->next != (*f_w))
	{
		b_cour = (*f_w)->next;
		(*f_w)->next = b_cour->next;
		delete b_cour;
	}

	delete (*f_w);

	(*f_w) = NIL_BN;
}




/*Print the set of b_Node W*/
void Print_b_Node_Set(b_Node *first_w,FILE *out_file)
{
	b_Node *n_cour;

	n_cour = first_w;
	while(n_cour != NIL_BN)
	{
		fprintf(out_file,"%ld ",n_cour->id);

		n_cour = n_cour->next;
	}
	fputs("",out_file);
}


void Print_LCC_b_Node(b_Node *first_w,FILE *out_file)
{
	b_Node *n_cour;

	n_cour = first_w;
	do
	{
		fprintf(out_file,"%ld ",n_cour->id);

		n_cour = n_cour->next;
	}
	while(n_cour->id != first_w->id);
	fputs("",out_file);
}



/*Print the b_Edge set delta(W)*/
void Print_b_Edge_Set(b_Edge *f_e,FILE *out_file)
{
	b_Edge *e_cour;

	e_cour = f_e;
	while(e_cour != NIL_BE)
	{
		fprintf(out_file,"%ld %ld   %f  E%ld\n",e_cour->node1,e_cour->node2,e_cour->cap,e_cour->num);

		e_cour = e_cour->next;
	}

	fputs("",out_file);
}

/*This function renumbers the nodes of a given set of nodes W, in order to apply*/
/*the ghct function (for instance) on the subgraph induced by W (ie G[W]). You must*/
/*notice that the ghct function and the NOI function work on graph which nodes are*/
/*consecutives (ie numerotate from 1 to n).*/
/*f_w is the head of the linked list of nodes representing the set W.*/
/*Tab1_2 gives the the correspondant of a node of W in the new numerotation.*/
/*n_12 is the number of nodes of the graph G from which the set W is comming from.*/
/*Tab_2_1 allows to make the correspondance in the back way.*/
/*n_21 is the number of nodes in the set W.*/
void RenumberTables(b_Node *f_w,simpleNode **Tab_12,long n_12,simpleNode **Tab_21,long *n_21)
{
	b_Node *b_cour;
	simpleNode *aux;
	long k;

	/*Allocation of Tab_12 and Tab_21:the default value will be 0*/
	(*Tab_12) = new simpleNode[n_12];
	(*Tab_21) = new simpleNode[n_12];

	b_cour = f_w;
	*n_21 = 0;
	while(b_cour != NIL_BN)
	{
		(*n_21)++;

		k = b_cour->id;
		(*Tab_12)[k-1] = (*n_21);

		(*Tab_21)[(*n_21)-1] = k;

		b_cour = b_cour->next;
	}

	aux = (*Tab_21);
	(*Tab_21) = new simpleNode[(*n_21)];

	for(k=0;k<(*n_21);k++)
	{
		(*Tab_21)[k] = aux[k];
	}

	delete[] aux;
}


/*Create a new list of b_Node from the list f_w. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
b_Node *renumber_nodes(b_Node *f_w,simpleNode *corresp_tab,long n_corresp)
{
	b_Node *b_cour,*f_new,**new_cour;
	long k;

	f_new = NIL_BN;
	new_cour = &f_new;

	b_cour = f_w;
	while(b_cour != NIL_BN)
	{
		k = b_cour->id;

		(*new_cour) = new b_Node;
		(*new_cour)->id = corresp_tab[k-1];
		(*new_cour)->next = NIL_BN;

		b_cour = b_cour->next;
		new_cour = &((*new_cour)->next);
	}

	return f_new;
}



/*Create a new list of simpleEdge from the list Tab. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
simpleEdge *renumber_edges(simpleEdge *Tab,long n_tab,simpleNode *corresp_tab,long n_corresp)
{
	long i,k;
	simpleEdge *tab_new;

	tab_new = new simpleEdge[n_tab];

	for(i=0;i<n_tab;i++)
	{
		k = Tab[i].node1;
		tab_new[i].node1 = corresp_tab[k-1];

		k = Tab[i].node2;
		tab_new[i].node2 = corresp_tab[k-1];

		tab_new[i].cap = Tab[i].cap;
	}

	return tab_new;
}




/*This function similar to renumber_nodes() and renumber_edges().*/
/*The new function achieve in one the same time the job done by the*/
/*renumber_nodes(), renumber_edges() and RenumberTables().*/
/*It returns the new list of edges obtained by renumbering.*/
simpleEdge *renumber_edge_list(simpleEdge *Te,long m,simpleNode **T1,long n1,simpleNode **T2,long *n2)
{
	long i,node1,node2,new_number1,new_number2;
	simpleEdge *Tab;


	Tab = new simpleEdge[m];

	(*T1) = new simpleNode[n1];
	(*T2) = new simpleNode[n1];


	(*n2) = 0;

	for(i=0;i<n1;i++)
	{
		(*T1)[i] = 0;
		(*T2)[i] = 0;
	}

	//cout << "M == " << m << endl;

	for(i=0;i<m;i++)
	{
		node1 = Te[i].node1;
		node2 = Te[i].node2;

		//cout << "Node1 == " << node1 << " " << "Node2 == " << node2 << endl;

		/*If the node has already been relabel then we take its new label.*/
		/*Else we insert a new label for the node in the new labels table T1*/
		/*and we insert the node in the old labels table T2.*/
		if((*T1)[node1-1] != 0)
		{
			new_number1 = (*T1)[node1-1];
		}
		else
		{
			(*n2)++;
			(*T1)[node1-1] = (*n2);
			(*T2)[(*n2)-1] = node1;

			new_number1 = (*n2);
		}


		/*Same as node1*/
		if((*T1)[node2-1] != 0)
		{
			new_number2 = (*T1)[node2-1];
		}
		else
		{
			(*n2)++;
			(*T1)[node2-1] = (*n2);
			(*T2)[(*n2)-1] = node2;

			new_number2 = (*n2);
		}


		/*We insert a new edge in the new edge table*/
		Tab[i].node1 = new_number1;
		Tab[i].node2 = new_number2;
		Tab[i].cap = Te[i].cap;
	}

	return Tab;
}





/*Inserts a b_Node in the given b_Node linked list.*/
/*The function inserts the b_Node such as the first*/
/*node of the list is the one with the lower id.*/
void insert_b_Node(b_Node **f_list,b_Node **b_cour)
{
	if(*f_list != NIL_BN)
	{
		if((*f_list)->id > (*b_cour)->id)
		{
			(*b_cour)->next = (*f_list);
			(*f_list) = (*b_cour);
		}
		else
		{
			(*b_cour)->next = (*f_list)->next;
			(*f_list)->next = (*b_cour);
		}
	}
	else
	{
		(*b_cour)->next = (*f_list);
		(*f_list) = (*b_cour);
	}
}





double max(double a, double b)
{
 	if (a>b) return a;
  	else return b;
}

double min(double a, double b)
{
  	if (a<b) return a;
  	else return b;
}


/*Renvoie l'indice dans la liste de l'arete e = (node1,node2)*/
long search_edge(simpleEdge *eList,long n_List,long node1,long node2)
{
	long i;
	Bool trouver;

	trouver = false; /// instead of defined False
	i = -1;
	do
	{
		i++;

		if( ((eList[i].node1 == node1) && (eList[i].node2 == node2)) ||
			((eList[i].node1 == node2) && (eList[i].node2 == node1)) )
		{
			trouver = true;  /// instead of defined True
		}
	}
	while((i < n_List-1) && (!trouver));

	if(trouver == true)  /// instead of defined True
	{
		return i;
	}
	else return -1;
}


/*Calcul l'intersection des cercles (x1,y1,R1) et (x2,y2,R2)*/
/*Renvoie:						    */
/* +   -1 si le calcul est impossible.*/
/* +   0 s'il n'y a pas d'intersection.*/
/* +   1 s'il y a un seul point.Dans ce cas le point est stocke dans (xa,ya).*/
/* +   2 s'il y a deux points d'intersection.*/
int intersection_cercle(double x1,double y1,double R1,double x2,double y2,double R2,double
*xa,double *ya,double *xb,double *yb)
{
	double cN,cA,cB,cC,delta;

	/*Impossible de calculer avec cette formule*/
	if(y1 == y2)
		return -1;

	cN = (R2*R2 - R1*R1 - x2*x2 + x1*x1 - y2*y2 + y1*y1)/(2*(y1 - y2));
	cA = pow((x1 - x2)/(y1 - y2),2) + 1.0;
	cB = 2*y1*(x1 - x2)/(y1 - y2) - 2*cN*(x1 - x2)/(y1 - y2) - 2*x1;
	cC = x1*x1 + y1*y1 + cN*cN - R1*R1 - 2*y1*cN;

	delta = cB*cB - 4*cA*cC;

	/*Pas d'intersection pour ces deux cercles*/
	if(delta < 0)
		return 0;

	/*Dans ce cas il y a un seul point d'intersection*/
	if(delta == 0)
	{
		*xa = (-cB + sqrt(delta))/(2*cA);
		*ya = cN - (*xa)*(x1 - x2)/(y1 - y2);

		return 1;
	}

	/*Dans ce cas il y a deux points d'intersection*/
	if(delta > 0)
	{
		*xa = (-cB + sqrt(delta))/(2*cA);
		*ya = cN - (*xa)*(x1 - x2)/(y1 - y2);

		*xb = (-cB - sqrt(delta))/(2*cA);
		*yb = cN - (*xb)*(x1 - x2)/(y1 - y2);

		return 2;
	}

	return 0;
}


/*Recherche dans une liste de b_Node le b_Node d'ID id*/
/*La fonction renvoie le précédent sur l'élément*/
b_Node **recherche_b_node(b_Node **list,long id)
{
	b_Node **prec,*cour;

	prec = list;
	cour = (*list);
	while((cour != NIL_BN) && (cour->id <= id))
	{
		prec = &(cour->next);
		cour = cour->next;
	}

	return prec;
}


//On s'assure que le sommet de plus petit indice
//soit mis en premier
void Insere_Node(b_Node **W,long id)
{
	b_Node *cour;

	cour = new b_Node;
	cour->id = id;
	cour->next = NIL_BN;

	if((*W) != NIL_BN)
	{
		if((*W)->id < cour->id)
		{
			cour->next = (*W)->next;
			(*W)->next = cour;
		}
		else
		{
			cour->next = (*W);
			(*W) = cour;
		}
	}
	else
	{
		(*W) = cour;
	}
}

//Ici, on insere les noeuds de sorte que
//s'il y a un terminal, il est mis en premier
//sinon, le plus petit sommet est mis en premier

void Insere_Node(b_Node **W,long id,bool terminal)
{
	b_Node *cour;

	cour = new b_Node;
	cour->id = id;
	cour->next = NIL_BN;

	if((*W) != NIL_BN)
	{
		if(terminal)
		{
			cour->next = (*W);
			(*W) = cour;
		}
		else
		{
			cour->next = (*W)->next;
			(*W)->next = cour;
		}
	}
	else
	{
		(*W) = cour;
	}
}


void Insere_Edge(b_Edge **L,b_Edge *be)
{
	if((*L) != NIL_BE)
	{
		be->next = (*L);
		(*L) = be;
	}
	else
	{
		(*L) = be;
	}
}
