#include <iostream>
#include "I_digraphe.h"

bool operator==(I_dinode &n1,I_dinode &n2)
{
	return (n1.num_reel == n2.num_reel);
}

bool operator!=(I_dinode &n1,I_dinode &n2)
{
	return (n1.num_reel != n2.num_reel);
}

bool operator<(I_dinode &n1,I_dinode &n2)
{
	return (n1.num_reel < n2.num_reel);
}

bool operator>(I_dinode &n1,I_dinode &n2)
{
	return (n1.num_reel > n2.num_reel);
}

ostream &operator<<(ostream &os,I_dinode &n1)
{
	os << n1.num_reel;
	return os;
}

ostream &operator<<(ostream &os,I_diarc &a)
{
	os << "(" << *a.source << "->" << *a.dest << ",";
	os << a.cap << "," << a.num << ")";
	return os;
}


I_digraph::I_digraph(int s,int t,I_Graph *g,int k)
{
	int ms,mt,muv;

	n = g->n_Nodes;
	m = 0;

	//Source et destination
	source.num = 1;
	source.num_reel = s;
	source.premier = (I_diarc *)0;

	dest.num = 2;
	dest.num_reel = t;
	dest.premier = (I_diarc *)0;

	//Construction de V1 et V2
	V1 = new I_dinode[n-2];
	V2 = new I_dinode[n-2];

	int i,j=0;
	for(i=1;i<=n;i++)
	{
		if((i != s) && (i != t))
		{
			V1[j].num = j+3;
			V1[j].num_reel = i;
			V1[j].premier = (I_diarc *)0;

			V2[j].num = j+3;
			V2[j].num_reel = i + n;		//les sommets de V2 sont identifies par i+n
			V2[j].premier = (I_diarc *)0;

			j++;
		}
	}

	if(g == (I_Graph *)0)
		return;

	//Calcul de degre de s
	ms = 0;
	I_Edge *ecour = g->Nodes[s-1].first_edge;
	while(ecour != NIL_E)
	{
		ms++;
		ecour = ecour->next;
	}


	//Calcul de degre de t
	mt = 0;
	ecour = g->Nodes[t-1].first_edge;
	while(ecour != NIL_E)
	{
		//S'il existe une arete entre s et t
		//il faut eviter de la compter 2 fois
		if(ecour->adjac->id != s)
			mt++;

		ecour = ecour->next;
	}

	//Le nombre des aretes	uv, u,v != s,t
	muv = g->m_Edges - ms - mt;

//	cout << "degre de s = " << ms << endl;
//	cout << "degre de t = " << mt << endl;
//	cout << "degre de uv = " << muv << endl;

	//Construction de la liste des arcs
	//On va parcourir chaque sommet du graphe non oriente
	//et inserer les arcs correspondant a chaque arete
	//de la liste d'adjacence du sommet

	//Allocation du tableau des arcs
	larc = new I_diarc[ms+mt+2*muv+k*(n-2)];

	int u,v;
	m = 0;
	I_dinode no;

	//Construction des arcs de la source
	ecour = g->Nodes[source.num_reel-1].first_edge;
	while(ecour != NIL_E)
	{
		larc[m].num = m+1;
		larc[m].source = &source;
		larc[m].f = ecour->X;
		larc[m].cap = 1;
		larc[m].suivant = (I_diarc *)0;
		larc[m].pedge = ecour;

		if(ecour->adjac->id != dest.num_reel)
		{
			no.num = -1;
			no.num_reel = ecour->adjac->id;
			no.premier = (I_diarc *)0;

			u = rechercheDichotomique(V1,no,0,n-3);
			larc[m].dest = &V1[u];
		}
		else
		{
			larc[m].dest = &dest;
		}

		lierArc(source,&larc[m]);

		ecour = ecour->next;
		m++;
	}

	//Construction des arcs en direction de f
	//Les sommets qui auront un arc vers t
	//sont en fait ceux qui sont dans la liste
	//d'adjacence de t dans le graphe oriente.
	//Sauf evidemment le sommet s.

	//Construction des arcs de la destination
	ecour = g->Nodes[dest.num_reel-1].first_edge;
	while(ecour != NIL_E)
	{
		if(ecour->adjac->id != source.num_reel)
		{
			no.num = -1;
			no.num_reel = ecour->adjac->id+n;
			no.premier = (I_diarc *)0;

			v = rechercheDichotomique(V2,no,0,n-3);
			larc[m].source = &V2[v];

			larc[m].num = m+1;
			larc[m].dest = &dest;
			larc[m].f = ecour->cap;
			larc[m].cap = 1;
			larc[m].suivant = (I_diarc *)0;
			larc[m].pedge = ecour;

			lierArc(V2[v],&larc[m]);

			m++;
		}

		ecour = ecour->next;
	}

	//Construction des arcs (u,v') et (v,u')
	//A partir des aretes uv, u,v != s,t

	for(i=1;i<=n;i++)
	{
		if((i != s) && (i != t))
		{
			ecour = g->Nodes[i-1].first_edge;

			no.num = -1;
			no.num_reel = i;
			no.premier = (I_diarc *)0;

			u = rechercheDichotomique(V1,no,0,n-3);

			//On ajoute les arcs correspondant a chaque voisin
			while(ecour != NIL_E)
			{
				if((ecour->adjac->id != s) && (ecour->adjac->id != t))
				{
					no.num = -1;
					no.num_reel = ecour->adjac->id+n;
					no.premier = (I_diarc *)0;

					v = rechercheDichotomique(V2,no,0,n-3);

					larc[m].num = m+1;
					larc[m].source = &V1[u];
					larc[m].dest = &V2[v];
					larc[m].f = ecour->cap;
					larc[m].cap = 1;
					larc[m].suivant = (I_diarc *)0;
					larc[m].pedge = ecour;

					lierArc(V1[u],&larc[m]);

					m++;
				}

				ecour = ecour->next;
			}
		}
	}

	//On rajoute les aretes paralleles (u,u')
	for(i=0;i<n-2;i++)
	{
		for(j=0;j<k;j++)
		{
			larc[m].num = m+1;
			larc[m].source = &V1[i];
			larc[m].dest = &V2[i];
			larc[m].f = 0;
			larc[m].cap = 1;
			larc[m].suivant = (I_diarc *)0;
			larc[m].pedge = NIL_E;

			lierArc(V1[i],&larc[m]);

			m++;
		}
	}
}

I_digraph::~I_digraph()
{
	delete[] V1;
	delete[] V2;
	delete[] larc;

	V1 = V2 = (I_dinode *)0;
	larc = (I_diarc *)0;
}


int I_digraph::getNbNode()
{
	return (2*n-2);
}

int I_digraph::getNbArc()
{
	return m;
}

bool I_digraph::lierArc(I_dinode &sommet,I_diarc *a)
{
	if(sommet.premier == (I_diarc *)0)
		sommet.premier = a;
	else
	{
		a->suivant = sommet.premier;
		sommet.premier = a;
	}

	return true;
}

/*bool I_digraph::addArc(I_dinode &sommet,I_diarc &a)
{
	larc.push_back(a);
	lierArc(sommet,&larc[larc.size()-1]);

	return true;
}*/


void I_digraph::printVoisin(ostream &os,I_dinode &n1)
{
	I_diarc *a;

	a = n1.premier;
	while(a != (I_diarc *)0)
	{
		os << *a << endl;

		a = a->suivant;
	}
}


void I_digraph::print(ostream &os)
{
	printVoisin(os,source);

	for(int i=0;i<n-2;i++)
	{
		printVoisin(os,V1[i]);
	}

	for(int i=0;i<n-2;i++)
	{
		printVoisin(os,V2[i]);
	}
}

void I_digraph::printTable(ostream &os)
{
	for(int i=0;i<m;i++)
	{
		os << larc[i] << endl;
	}
}


void I_digraph::graphViz(const char *nomfich,bool avec_val)
{
	ofstream fich(nomfich);

	fich << "I_digraph G{" << endl;

	/*Fixation de la position des sommets du label*/
	fich << "node [shape=circle];" << endl;

	I_Node *nptr,*n1;
	char s[256];
	string sh_node_label = "";

	fich << source << " [label = \"" << source << "\"]" << endl;
	fich << dest << " [label = \"" << dest << "\"]" << endl;

	for(int i=0;i<n-2;i++)
	{
		fich << V1[i] << " [label = \"" << V1[i] << "\"]" << endl;
		fich << V2[i] << " [label = \"" << V2[i].num_reel-n << "\"]" << endl;
	}


	I_diarc *a = source.premier;

	while(a != NIL_A)
	{
		if(a->cap >= 1.0)
		{
			fich << *a->source << "->" << *a->dest << "[id=\"";
			fich << a->num << "\"];" << endl;
		}
		else
		{
			if(a->cap > 0.0)
			{
				if(avec_val == true)
				{
					fich << *a->source << "->" << *a->dest << "[id=\"" << a->num << "\"]";
					fich << "[style=\"dashed\"][color=\"red\"]";
					fich << "[label=\"" << a->cap << "\"];" << endl;
				}
				else
				{
					fich << *a->source << "->" << *a->dest;
					fich << "[id=\"" << a->num << "\"]";
					fich << "[style=\"dashed\"][color=\"red\"];" << endl;
				}
			}
		}

		a = a->suivant;
	}


	for(int i=0;i<n-2;i++)
	{
		a = V1[i].premier;

		while(a != NIL_A)
		{
			if(a->cap >= 1.0)
				fich << *a->source << "->" << *a->dest << "[id=\"" << a->num << "\"];" << endl;
			else
			{
				if(a->cap > 0.0)
				{
					if(avec_val == true)
					{
						fich << *a->source << "->" << *a->dest << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"]";
						fich << "[label=\"" << a->cap << "\"];" << endl;
					}
					else
					{
						fich << *a->source << "->" << *a->dest;
 						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"];" << endl;
					}
				}
			}

			a = a->suivant;
		}
	}

	for(int i=0;i<n-2;i++)
	{
		a = V2[i].premier;

		while(a != NIL_A)
		{
			if(a->cap >= 1.0)
				fich << *a->source << "->" << *a->dest << "[id=\"" << a->num << "\"];" << endl;
			else
			{
				if(a->cap > 0.0)
				{
					if(avec_val == true)
					{
						fich << *a->source << "->" << *a->dest << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"]";
						fich << "[label=\"" << a->cap << "\"];" << endl;
					}
					else
					{
						fich << *a->source << "->" << *a->dest;
 						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"];" << endl;
					}
				}
			}

			a = a->suivant;
		}
	}



	fich << "}" << endl;

	fich.close();




}


void I_digraph::graphViz(const char *nomfich,simpleEdge *solution,bool avec_val)
{
	ofstream fich(nomfich);

	fich << "I_digraph G{" << endl;

	/*Fixation de la position des sommets du label*/
	fich << "node [shape=circle];" << endl;

	char s[256];
	string sh_node_label = "";

	fich << source << " [label = \"" << source << "\"]" << endl;
	fich << dest << " [label = \"" << dest << "\"]" << endl;

	for(int i=0;i<n-2;i++)
	{
		fich << V1[i] << " [label = \"" << V1[i] << "\"]" << endl;
		fich << V2[i] << " [label = \"" << V2[i].num_reel-n << "'" << "\"]" << endl;
	}

	I_diarc *a = source.premier;

	while(a != NIL_A)
	{
		if(solution[a->num-1].cap >= 1.0)
		{
			fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2 << "[id=\"";
			fich << a->num << "\"];" << endl;
		}
		else
		{
			if(solution[a->num-1].cap > 0.0)
			{
				if(avec_val == true)
				{
					fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
					fich << "[id=\"" << a->num << "\"]";
					fich << "[style=\"dashed\"][color=\"red\"]";
					fich << "[label=\"" << solution[a->num-1].cap << "\"];" << endl;
				}
				else
				{
					fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
					fich << "[id=\"" << a->num << "\"]";
					fich << "[style=\"dashed\"][color=\"red\"];" << endl;
				}
			}
		}

		a = a->suivant;
	}


	for(int i=0;i<n-2;i++)
	{
		a = V1[i].premier;

		while(a != NIL_A)
		{
			if(solution[a->num-1].cap >= 1.0)
			{
				fich << solution[a->num-1].node1 << "->";
				fich << solution[a->num-1].node2;
				fich << "[id=\"" << a->num << "\"];" << endl;
			}
			else
			{
				if(solution[a->num-1].cap > 0.0)
				{
					if(avec_val == true)
					{
						fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"]";
						fich << "[label=\"" << a->cap << "\"];" << endl;
					}
					else
					{
						fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
 						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"];" << endl;
					}
				}
			}

			a = a->suivant;
		}
	}

	for(int i=0;i<n-2;i++)
	{
		a = V2[i].premier;

		while(a != NIL_A)
		{
			if(solution[a->num-1].cap >= 1.0)
			{
				fich << solution[a->num-1].node1 << "->";
				fich << solution[a->num-1].node2;
				fich << "[id=\"" << a->num << "\"];" << endl;
			}
			else
			{
				if(solution[a->num-1].cap > 0.0)
				{
					if(avec_val == true)
					{
						fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"]";
						fich << "[label=\"" << a->cap << "\"];" << endl;
					}
					else
					{
						fich << solution[a->num-1].node1 << "->" << solution[a->num-1].node2;
 						fich << "[id=\"" << a->num << "\"]";
						fich << "[style=\"dashed\"][color=\"red\"];" << endl;
					}
				}
			}

			a = a->suivant;
		}
	}



	fich << "}" << endl;

	fich.close();
}


void I_digraph::getDelta(list<I_dinode *> &W,list<I_diarc *> &delta_w)
{
	bool *marquer;
	marquer = new bool[2*n];

	list<I_dinode *>::iterator it;

	for(int i=0;i<2*n;i++)
	{
		marquer[i] = false;
	}

	for(it=W.begin();it!=W.end();it++)
	{
		marquer[(*it)->num_reel-1] = true;
	}

	for(it=W.begin();it!=W.end();it++)
	{
		I_diarc *a;

		a = (*it)->premier;

		//cout << " -- " << *(*it) << endl;

		while(a != NIL_A)
		{
			//cout << " - " << *a->dest << endl;

			if(!marquer[a->dest->num_reel-1])
			{
				//cout << " - " << *a->dest << endl;
				delta_w.push_back(a);
			}

			a = a->suivant;
		}
	}

	delete[] marquer;
}

/*void I_digraph::getDeltaPlusS(int s,list<I_diarc *> &delta_s,int edge_flag)
{
	int u = findSommet(source,s);

	I_diarc *a = source[u]->premier;

	while(a != NIL_A)
	{
		switch(edge_flag)
		{
			case ALL_ARCS: 	delta_s.push_back(a);
					break;

			case NO_0_ARCS:	if(a->cap > 0)
						delta_s.push_back(a);
					break;

			default:	cout << "Delta+S() choix par defaut" << endl;
					delta_s.push_back(a);
					break;
		}

		a = a->suivant;
	}
}

void I_digraph::getDeltaMoinsT(int t,list<I_diarc *> &delta_t,int edge_flag)
{
	I_diarc *a;

	for(int i=0;i<n;i++)
	{
		a = V2[i]->premier;

		while(a != NIL_A)
		{
			if(a->dest->num_reel == (t+3*n))
			{
				switch(edge_flag)
				{
					case ALL_ARCS: 	delta_t.push_back(a);
							break;

					case NO_0_ARCS:	if(a->cap > 0)
								delta_t.push_back(a);
							break;

					default:	cout << "Erreur : Delta-T() choix par defaut" << endl;
							delta_t.push_back(a);
							break;
				}

			}

			a = a->suivant;
		}
	}

	//Il ne faut pas oublier d'ajouter un eventuel arc (t',t)
	a = V1[t-1]->premier;
	while(a != NIL_A && a->dest->num_reel != t+3*n)
		a = a->suivant;

	if(a != NIL_A)
	{
		switch(edge_flag)
		{
			case ALL_ARCS: 	delta_t.push_back(a);
					break;

			case NO_0_ARCS:	if(a->cap > 0)
						delta_t.push_back(a);
					break;

			default:	cout << "Erreur : Delta-T() choix par defaut" << endl;
					delta_t.push_back(a);
					break;
		}
	}
}*/


double I_digraph::makeFeasible(I_Graph *g,int k)
{
    //Ici, on va rendre realisable la solution courante pour la demande (s,t)
    //Idee:
    //  1-  Pour chaque demande (s,t)
    //      a-  Parcourir les X et mettre 1 sur tous les arcs dont le X vaut 1
    //          et mettre EPS sur tous les arcs dont le X est fractionnaire (pour pouvoir les retrouver plus facilement)
    //
    //      b-  calculer un flot max pour déterminer le nombre de chemins arc-disjoints
    //          qu'on obtient avec les arcs à 1
    //      c-  mettre le flot à 0 sur les arcs à EPS, leur donner une capacité 1
    //          et mettre à 1 le flot sur tous les arcs entiers
    //      d-  Tantque Rst < k Faire
    //              trouver une chaine augmentantes entre s et t
    //              mettre à jour le chemin dans le graphe
    //              mettre à jour Rst
    //          Fin
    //      e-  Pour tous les arcs du graphe qui ont un flot > 0 marquer dans G le X correspondant à 1

    double MyEPS = 0.0001;

    //1ere MAJ des capacités du graphe
    for(int j=0;j<getNbArc();j++)
    {
        if((larc[j].pedge != NIL_E) && (larc[j].pedge->X > 0))
        {
            if(larc[j].pedge->X == 1)
            {
                    larc[j].cap = 1;
            }
            else
            {
                larc[j].cap = MyEPS;
            }
        }
    }

    //Calcul du nombre de chemin arc-disjoint dejà présents
    int s_p,t_p;

    bool *marquer = new bool[2*n];

    I_Graph_Flot *gflot = new I_Graph_Flot;
    directInitGraphFlot(gflot,2*n);

    s_p = source.num_reel-1;
    t_p = dest.num_reel-1;

    double flot_max = gflot->flotMaxAugmentation(s_p,t_p,marquer);
    //cout << "Nb Chemin = " << flot_max << endl;

    int nbChemin = (int)flot_max;

    //2e MAJ des capacite des arcs du graphe
    I_Edge_Flot *eflot;
    for(int j=0;j<gflot->size;j++)
    {
        eflot = gflot->A[j];
        while(eflot != ((I_Edge_Flot *) 0))
        {
            if(eflot->c > 0)
            {
                //Comme on ramène à 0 le flot sur les arcs EPS
                //il faut mettre à 1 le flot sur les arcs qui
                //forment un chemin avec des arcs EPS
                //Attention: certains arcs entiers peuvent être
                //utilisés seulement avec des arcs EPS.
                //Pour différencier les vrais arcs entiers à 1
                //et les faux, on regarde le flot.
                //Comme EPS est suffisament petit, si le flot est > 0.5
                //alors c'est un vrai arc à 1. Donc on met son flot à 0
                //si le flot est < 0.5 alors c'est un faux arc à 1. Donc on met
                //son flot à 0.

                //cout << "eflot->c = " << eflot->c << " My EPS = " << MyEPS << endl;

                if(eflot->c <= MyEPS)
                {
                    eflot->f = 0;
                    eflot->c = 1;
                }
                else if(eflot->f > 0.5)
                    eflot->f = 1;
                else
                    eflot->f = 0;
            }

            eflot = eflot->next;
        }
    }

    //Calcul des nouvelles chaines augmentantes
    double augm;
	//cout << "NbChem == " << nbChemin << endl;

	do
	{
		//Calcul une chaine augmentante et met à jour le flot sur les arcs conrcernés
		augm = gflot->AugmenteFlot(s_p,t_p,marquer);

		//cout << "Aug = " << augm << endl;

		if(augm == 1)
			nbChemin++;

	}while(augm != 0 && nbChemin < k);

    /*double augm;
    while(nbChemin < k)
    {
        //Calcul une chaine augmentante et met à jour le flot sur les arcs conrcernés
        augm = gflot->AugmenteFlot(s_p,t_p,marquer);

        //cout << "Augm = " << augm << endl;

        if(augm > 0)
            nbChemin++;
    }*/

    //3e MAJ du graphe
    tranfertFlot(gflot);

    for(int j=0;j<getNbArc();j++)
    {
        if(larc[j].f > 0)
        {
            larc[j].cap = 1;

            if(larc[j].pedge != NIL_E)
            {
                larc[j].pedge->X = 1;
                larc[j].pedge->back->X = 1;
            }
        }
    }

    //Suppression du graphe flot
    delete[] marquer;
    delete gflot;
}


void I_digraph::directInitGraphFlot(I_Graph_Flot *gr,int nodes)
{
  	int u, v;
  	double cap_inf,cap_sup;

	gr->InitGraph(nodes);

    I_diarc *a = source.premier;
    while(a != NIL_A)
    {
        a->f = 0;

        if(a->cap > 0)
        {
            u = a->source->num_reel;
            v = a->dest->num_reel;
            cap_inf = 0;
            cap_sup = a->cap;

            gr->AddEdge(u-1,v-1,cap_inf,cap_sup,a->num);
        }

        a = a->suivant;
    }

	for(int i=0;i<n-2/*V1.size()*/;i++)
	{
		I_diarc *a = V1[i].premier;
		while(a != NIL_A)
		{
			a->f = 0;

			if(a->cap > 0)
			{
				u = a->source->num_reel;
    			v = a->dest->num_reel;
    			cap_inf = 0;
   				cap_sup = a->cap;

    			gr->AddEdge(u-1,v-1,cap_inf,cap_sup,a->num);
			}

			a = a->suivant;
		}
	}

	for(int i=0;i<n-2/*V2.size()*/;i++)
	{
		I_diarc *a = V2[i].premier;
		while(a != NIL_A)
		{
			a->f = 0;

			if(a->cap > 0)
			{
				u = a->source->num_reel;
   				v = a->dest->num_reel;
   				cap_inf = 0;
   				cap_sup = a->cap;

   				gr->AddEdge(u-1,v-1,cap_inf,cap_sup,a->num);
			}

			a = a->suivant;
		}
	}

  	for (int i=0;i<nodes;i++)
   		gr->AddVertex(i);
}

//Transfert du flot graphe Graph_Flot vers le I_digraph
void I_digraph::tranfertFlot(I_Graph_Flot *gr)
{
	for(int i=0;i<gr->size;i++)
	{
		I_Edge_Flot *ecour;
		ecour = gr->A[i];

		while(ecour != (I_Edge_Flot *)0)
		{
			if(ecour->c > 0)
			{
				int num = ecour->numDiarc;
				if(ecour->f > EPS)
					larc[num-1].f = ecour->f;
				else
					larc[num-1].f = 0;

				//larc[num-1]->label = ecour->label;
			}

			ecour = ecour->next;
		}
	}
}


double I_digraph::calculFlot(list<I_dinode *> &W)
{
	int s_p,t_p;

	s_p = source.num_reel-1;
	t_p = dest.num_reel-1;

	bool *marquer = new bool[2*n];
	I_Graph_Flot *gflot = new I_Graph_Flot;


	directInitGraphFlot(gflot,2*n);

	double flot_max = gflot->flotMaxAugmentation(s_p,t_p,marquer);
	tranfertFlot(gflot);

	/*gflot->GraphViz("sortie.dot");
	system("dot -Tps -o sortie.ps sortie.dot");*/

	//Construction de l'ensemble W

    //La source
    if(marquer[s_p])
    {
        W.push_back(&source);
    }

	//On regarde d'abord si la premiere courche sont marquees
	for(int i=0;i<n-2;i++)
	{
		if(marquer[V1[i].num_reel-1])
			W.push_back(&V1[i]);
	}

	for(int i=0;i<n-2;i++)
	{
		if(marquer[V2[i].num_reel-1])
			W.push_back(&V2[i]);
	}

	delete[] marquer;

	//Suppression du graphe flot
	delete gflot;

	return flot_max;
}


bool corresponde(I_diarc *a1,I_diarc *a2,int n)
{
    /*cout << "A1 = " << *a1 << endl;
    cout << "A2 = " << *a2 << endl;*/

	int p1,q1,p2,q2;

	if(a1 != NIL_A && a2 != NIL_A)
	{
		p1 = a1->source->num_reel;
		q1 = a1->dest->num_reel;

		p2 = a2->source->num_reel;
		q2 = a2->dest->num_reel;

		/*cout << "P1 == " << p1 << " " << q1 << endl;
		cout << "P2 == " << p2 << " " << q2 << endl;
		cout << "xxxxxxx" << endl;*/

		if(p1 <= n)
			p1 = p1;
		else if(p1 <= 2*n)
			p1 = p1-n;

		if(q1 <= n)
			q1 = q1;
		else if(q1 <= 2*n)
			q1 = q1-n;

		if(p2 <= n)
			p2 = p2;
		else if(p2 <= 2*n)
			p2 = p2-n;

		if(q2 <= n)
			q2 = q2;
		else if(q2 <= 2*n)
			q2 = q2-n;

		/*cout << "P1 == " << p1 << " " << q1 << endl;
		cout << "P2 == " << p2 << " " << q2 << endl;*/

		return ((p1 == p2 && q1 == q2) || (p1 == q2 && q1 == p2));
	}
	else
	{
		cout << "Erreur dans Corresponde 2:" << "A1 ou A2 NULL" << endl;
		return false;
	}
}

void I_digraph::getDeltaPlusS(list<I_diarc *> &delta_s,int edge_flag)
{
	I_diarc *a = source.premier;

	while(a != NIL_A)
	{
		switch(edge_flag)
		{
			case ALL_ARCS: 	delta_s.push_back(a);
                            break;

			case NO_0_ARCS:	if(a->cap > 0)
                            delta_s.push_back(a);
                            break;

			default:	    cout << "Delta+S() choix par defaut" << endl;
                            delta_s.push_back(a);
                            break;
		}

		a = a->suivant;
	}
}

void I_digraph::getDeltaMoinsT(list<I_diarc *> &delta_t,int edge_flag)
{
    int t = dest.num_reel;

	I_diarc *a;

	for(int i=0;i<n-2;i++)
	{
		a = V2[i].premier;

		while(a != NIL_A)
		{
			if(a->dest->num_reel == t)
			{
				switch(edge_flag)
				{
					case ALL_ARCS: 	delta_t.push_back(a);
                                    break;

					case NO_0_ARCS:	if(a->cap > 0)
                                    delta_t.push_back(a);
                                    break;

                    default:	    cout << "Erreur : Delta-T() choix par defaut" << endl;
                                    delta_t.push_back(a);
                                    break;
				}

			}

			a = a->suivant;
		}
	}

	//Il ne faut pas oublier d'ajouter un eventuel arc (s,t)
    a = source.premier;
	while(a != NIL_A && a->dest->num_reel != t)
		a = a->suivant;

	if(a != NIL_A)
	{
		switch(edge_flag)
		{
			case ALL_ARCS: 	delta_t.push_back(a);
                            break;

			case NO_0_ARCS:	if(a->cap > 0)
                            delta_t.push_back(a);
                            break;

            default:	    cout << "Erreur : Delta-T() choix par defaut" << endl;
                            delta_t.push_back(a);
                            break;
		}
	}
}


bool I_digraph::arcLier2(I_diarc *a1,I_diarc *a2,int k,list<I_dinode *> &coupeW,list<I_diarc *> &lEntier)
{
    if(fabs(a1->cap+a2->cap-1) > EPS)
		return false;

    int s,t;

    s = source.num_reel;
    t = dest.num_reel;

    //Principe:
    //On met les arcs a1 et a2 à 0
    //on met une capacité infinie sur les arcs qui sont adjacents à a1 et a2
    //Puis on calcule un flot max entre s et t.

    //On modifie les capacites des arcs
    list<I_diarc *> lArcModif;
    list<double> lCapModif;

    lArcModif.push_back(a1);
    lCapModif.push_back(a1->cap);

    lArcModif.push_back(a2);
    lCapModif.push_back(a2->cap);

    a1->cap = 0;
    a2->cap = 0;

    //On met ensuite une capacité infinie sur les arcs adjacents
    //qui sortent de s et entrent sur a1 et a2

    //On commence par mettre +infini comme capacité sur les arcs fractionnaires de \delta+(s) et \delta-(t)
    I_diarc *a;
    list<I_diarc *> delta_s;
    getDeltaPlusS(delta_s,NO_0_ARCS);

    for(list<I_diarc *>::iterator it=delta_s.begin();it!=delta_s.end();it++)
    {
        a = *it;

        if(a->num != a1->num && a->num != a2->num && a->cap < 1)
        {
            lArcModif.push_back(a);
            lCapModif.push_back(a->cap);

            a->cap = INFINI;
        }
    }

    //On fait la même chose pour t
    list<I_diarc *> delta_t;
    getDeltaMoinsT(delta_t,NO_0_ARCS);

    for(list<I_diarc *>::iterator it=delta_t.begin();it!=delta_t.end();it++)
    {
        a = *it;

        if(a->num != a1->num && a->num != a2->num && a->cap < 1)
        {
            lArcModif.push_back(a);
            lCapModif.push_back(a->cap);

            a->cap = INFINI;
        }
    }

    //Si a1 est de la forme (u',v"), alors on interdit aussi
    //tous les arcs à 1 qui sont adjacents à a1. Même chose pour a2.
    //Evidemment, on ne le fait que si l'arc est à 1. Car sinon, il a déjà
    //été modifié.
    if((a1->source->num_reel != s) && (a1->source->num_reel <= n))
    {
        for(list<I_diarc *>::iterator it=delta_s.begin();it!=delta_s.end();it++)
        {
            a = *it;

            if(a->dest->num_reel == a1->source->num_reel && a->cap == 1)
            {
                lArcModif.push_back(a);
                lCapModif.push_back(a->cap);

                a->cap = INFINI;
            }
        }

        for(list<I_diarc *>::iterator it=delta_t.begin();it!=delta_t.end();it++)
        {
            a = *it;

            if(a->source->num_reel == a1->dest->num_reel && a->cap == 1)
            {
                lArcModif.push_back(a);
                lCapModif.push_back(a->cap);

                a->cap = INFINI;
            }
        }
    }

    if((a2->source->num_reel != s) && (a2->source->num_reel <= n))
    {
        for(list<I_diarc *>::iterator it=delta_s.begin();it!=delta_s.end();it++)
        {
            a = *it;

            if(a->dest->num_reel == a2->source->num_reel && a->cap == 1)
            {
                lArcModif.push_back(a);
                lCapModif.push_back(a->cap);

                a->cap = INFINI;
            }
        }

        for(list<I_diarc *>::iterator it=delta_t.begin();it!=delta_t.end();it++)
        {
            a = *it;

            if(a->source->num_reel == a2->dest->num_reel && a->cap == 1)
            {
                lArcModif.push_back(a);
                lCapModif.push_back(a->cap);

                a->cap = INFINI;
            }
        }
    }

    //Puis on calcule un flot max entre s et t
    coupeW.clear();
    double flotmax = calculFlot(coupeW);

    /*cout << "Coupe W = " << flotmax << endl;
    for(list<I_dinode *>::iterator it=coupeW.begin();it!=coupeW.end();it++)
    {
        //printNode(cout,*(*it));
        //cout << " ";
        cout << *(*it) << " ";
    }
    cout << endl;*/

    list<I_diarc *> delta_w;
    getDelta(coupeW,delta_w);

    //cout << "Delta W" << endl;
    double poids = 0;
    int nbEntier = 0;
    for(list<I_diarc *>::iterator it=delta_w.begin();it!=delta_w.end();it++)
    {
        //cout << *(*it) << endl;
        poids += (*it)->cap;

        if((*it)->cap == 1)
        {
            nbEntier++;
            lEntier.push_back(*it);
            //cout << "Entier = " << *(*it) << endl;
        }
    }

    /*cout << "Poids de la coupe = " << poids << endl;
    cout << "Nb Entier = " << nbEntier << endl;*/

    //On restore les valeurs des arcs modifiés avant de sortir
    list<double>::iterator itCap = lCapModif.begin();
    for(list<I_diarc *>::iterator it=lArcModif.begin();it!=lArcModif.end();it++)
    {
        a = *it;
        a->cap = *itCap;

        itCap++;
    }

    if(nbEntier == k-1 && poids == k-1)
    {
          return true;
    }
    else
    {
        //cout << "Coupe Non" << endl;
        return false;
    }
}

