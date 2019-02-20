#include "extra_cuts.h"

#include <omp.h>
#include <mutex>

//////////////////////////////////PRINTINGS////////////////////////////////////
//#define printings
//#define output_verif_cons
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////OUTPUT MPI////////////////////////////////////
#ifdef output_verif_cons
    string procOut_verif_cons_name;
    ofstream procOut_verif_cons;
#endif // output_verif_cons
///////////////////////////////////////////////////////////////////////////////


mutex l1,l2,l3;


extra_cuts::extra_cuts()
{
    //ctor
}


extra_cuts::~extra_cuts()
{
	delete[] x_coord;
	delete[] y_coord;

	//Destruction des graphes de demandes
	for(int i=0;i < vectDemande.size();i++)
	{
		delete grapheDemande[i];
	}

	delete[] grapheDemande;

	//Suppression du graphe non oriente
	DeleteGraph(&G);
	DeleteGraph(&G_dem);
}


simpleEdge *extra_cuts::readEUC2D(ifstream &fich)
{
	int i,j;
	double dx,dy;
	double cap;
	simpleEdge *edge_list;
	char buffer[256];


	/*Lecture des coordonn�es de chaque points*/
	x_coord = new double[tsp_n_nodes];
	y_coord = new double[tsp_n_nodes];

	tsp_m_edges = (tsp_n_nodes*(tsp_n_nodes-1))/2;


	for(int k=0;k<tsp_n_nodes;k++)
	{
		fich >> i;
		fich >> x_coord[i-1];
		fich >> y_coord[i-1];
	}

	/*construction de la liste des ar�tes avec les capacit�s*/
	edge_list = new simpleEdge[tsp_m_edges];

	long k = 0;
	for(i=0;i<tsp_n_nodes-1;i++)
	{
		for(j=i+1;j<tsp_n_nodes;j++)
		{
			dx = x_coord[i] - x_coord[j];
			dy = y_coord[i] - y_coord[j];

			cap = floor(sqrt(dx*dx + dy*dy)+0.5);

			edge_list[k].node1 = i+1;
			edge_list[k].node2 = j+1;
			edge_list[k].cap = cap;

			k++;
		}
	}

	//Lecture des demandes
	do
	{
		fich >> buffer;

	}while((!fich.eof()) && (strcmp(buffer,"DEMAND_SECTION") != 0));

	int nb_dem,nSource,nDest,u;

	fich >> nSource;
	fich >> nDest;

	cout << "NS = " << nSource << " " << nDest << endl;

	for(int i=0;i<nSource;i++)
	{
		fich >> u;
		source.push_back(u);
	}

	for(int i=0;i<nDest;i++)
	{
		fich >> u;
		dest.push_back(u);
	}

	fich >> nb_dem;

	for(int i=0;i<nb_dem;i++)
	{
		demande d;

		fich >> d.source;
		fich >> d.dest;
		fich >> d.val;

		vectDemande.push_back(d);
	}

	return edge_list;
}



simpleEdge *extra_cuts::readTsp(char *file_name,int type)
{
	cout << "Fich " << file_name << endl;

	ifstream fich(file_name);
	char buffer[256];
	long i,j;
	bool typetrouve,dimensiontrouvee;
	simpleEdge *e_list;

	if(!fich)
	{
		cout << "Erreur. Le fichier que vous avez sp�cifiez n'existe pas" << endl;
//		exit(Fatal);
	}


	/*Lecture de la dimension*/
	do
	{
		fich >> buffer;

	}while(!fich.eof() && ( (strncmp(buffer,"DIMENSION",strlen("DIMENSION")) != 0) &&
		(strncmp(buffer,"DIMENSION:",strlen("DIMENSION:")) != 0)) );


	if(!fich.eof())
	{
		if( strcmp(buffer,"DIMENSION"/*,strlen("DIMENSION")*/) == 0 )
			fich >> buffer;

		dimensiontrouvee = true;
		fich >> tsp_n_nodes;
	}
	else
	{
		cout << "Erreur dans le fichier: la dimension n'est pas sp�cifi�e" << endl;
//		exit(Fatal);
	}


	/*Lecture du Type*/
	do
	{
		fich >> buffer;
	}while((!fich.eof()) && (strncmp(buffer,"EDGE_WEIGHT_TYPE",strlen("EDGE_WEIGHT_TYPE")) != 0) );


	if(!fich.eof())
	{
		/*fich >> buffer;
		if(strncmp(buffer,"EUC_2D",strlen("EUC_2D")) == 0)
			typetrouve = true;
		else
			typetrouve = false;*/
		typetrouve = true;
	}
	else
	{
		cout << "Erreur dans le fichier: le type euclidien n'est pas sp�cifi�e" << endl;
//		exit(Fatal);
	}

	if(!typetrouve)
	{
		cout << "Erreur. Le type n'a pas �t� trouv�" << endl;
//		exit(Fatal);
	}


	/*Lecture du Flag des coordonn�es*/
	if( (type != EXPLICIT_LOWER_DIAG_ROW) && (type != EXPLICIT_UPPER_DIAG_ROW) && (type !=
	EXPLICIT_LOWER_ROW) && (type != EXPLICIT_UPPER_ROW) )
	{
		do
		{
			fich >> buffer;
		}
		while(!fich.eof() && (strncmp(buffer,"NODE_COORD_SECTION",strlen("NODE_COORD_SECTION"))	!= 0) );
	}
	else
	{
		do
		{
			fich >> buffer;
		}
		while(!fich.eof() && (strncmp(buffer,"EDGE_WEIGHT_SECTION",strlen("EDGE_WEIGHT_SECTION")) != 0) );
	}

	if(fich.eof())
	{
		cout << "Erreur. Les flags NODE_COORD_SECTION et EDGE_WEIGHT_SECTION n'ont pas �t� trouv�s" << endl;
//		exit(Fatal);
	}

	if(type == EUC_2D)
		e_list = readEUC2D(fich);
//	else if(type == GEO)
//		e_list = readGeo(fich);
//	else if(type == EXPLICIT_LOWER_ROW)
//		e_list = readExplicitLowerRow(fich);
//	else if(type == EXPLICIT_LOWER_DIAG_ROW)
//		e_list = readExplicitLowerDiagRow(fich);
//	else if(type == EXPLICIT_UPPER_ROW)
//		e_list = readExplicitUpperRow(fich);
//	else if(type == EXPLICIT_UPPER_DIAG_ROW)
//		e_list = readExplicitUpperDiagRow(fich);
//	else if(type == ATT)
//		e_list = readATT(fich);


	fich.close();

	return e_list;
}



/***********************************************************/
/*Lance la s�paration de toutes les familles de contraintes valides*/
int extra_cuts::separation(deque < int * > &seperated_cons, int type, FILE *sortie)
{
	I_Graph *gr;

	gr = &G;

	n = gr->n_Nodes;
	m = gr->m_Edges;

//	int nbctepartpool=0,nbctesppartpool=0,nbtotalpool=0,nbctedcpool=0,nbctetpcpool=0;

	//cout << "Separate " << numero_noeud << " sur " << NbSub << " au total." << endl;

	// partie à réintégrer si on a une pool de contraintes par défaut
//	nbctepartpool =	constraintPoolSeparation(0,getMaster()->getPoolPartition());
//	nbctesppartpool = constraintPoolSeparation(0,getMaster()->getPoolSPPartition());
//	nbctedcpool = constraintPoolSeparation(0,getMaster()->getPoolDoubleCut());
//	nbctetpcpool = constraintPoolSeparation(0,getMaster()->getPoolTriplePathCut());
//
//	nbtotalpool = nbctepartpool + nbctesppartpool + nbctedcpool + nbctetpcpool;
//
//
//	if(nbtotalpool != 0)
//		return nbtotalpool;


	/***********************************************************************************/
	/*****  Si toutes les coupes sont v�rifi�es et que le point est fractionnaire  *****/
	/*****  alors on r�duit le graphe et on s�pare les nouvelles contraintes:      *****/
	/*****			partitions,F-Partitions,roues paires...                *****/
	/************************ S�paration nouvelles contraintes *************************/
	/***********************************************************************************/


	/************************************************************************/
	/*****************  R�duction du graphe:op�ration Theta  ****************/

	b_Node *v_p,**bn_cour,*b_cour1;
	long nombre_contraction,n_contract;
	b_Node *contract_tab[N_SUR_2];
	I_Node *nptr;
	int nbContraction;

	for(int i=0;i<(n/2)+1;i++)
	{
		contract_tab[i] = NIL_BN;
	}

	//cout << "R�duction du graphe" << endl;

	/********************/
	/*Separation normale*/

	//cout << "Separation normale" << endl;
	int resultat_separation = 0;

	if(resultat_separation == 0)
	{
		resultat_separation = separation_globale(seperated_cons, 1, stdout);
	}

//	cout << "Id: " << id() << " LP == " << lp()->value() << endl;

	//A la racine de l'arbre on enregistre la solution fractionnaire
//	static bool sol_avant_branchement_trouvee = false;
//	if((sol_avant_branchement_trouvee == false) && (resultat_separation == 0) && (id() == 1))
//	{
//		cout << "Branchement" << endl;
//
//		sol_avant_branchement_trouvee = true;
//
//		getMaster()->solAvantBranchement = lp()->value();
//
//		string s = getMaster()->getFileName();
//		s += ".avbr";
//
//		getMaster()->majBestSol();
//		getMaster()->printBestSol((char *)s.c_str());
//		getMaster()->printStat();
//
//		//exit(Fatal);
//	}

	//cout << "NbNoeud == " << numero_noeud << endl;
	//create_branch_file(Te,m);


	/*D�contraction total du graphe:on fait simplement*/
	/*l'initialisation de tous les champs*/
	for(int i=0;i<n-1;i++)
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

	save_contraction_info(gr,1);
	save_contraction_info(gr,2);

#ifdef printings
	cout << "\n\nResultat == " << resultat_separation << endl << endl;
#endif // printings

	return (resultat_separation);
//    return (0);
}



void extra_cuts::init(const char *filename,bool isTSPFile,int type,int k_val,bool demandeConnectee)
{
	char dir[256];

	if( (strlen(filename) > 1) && (filename[0] == '.') && (filename[1] == '/') )
	{
		sprintf(dir,"%s",getenv("PWD"));

		sprintf(File_Name,"%s/%s",dir,filename+2);
	}
	else
	{
		if(filename[0] == '/')
		{
			sprintf(File_Name,"%s",filename);
		}
		else
		{
			sprintf(dir,"%s",getenv("PWD"));
			sprintf(File_Name,"%s/%s",dir,filename);
		}
	}

//	strncpy(OutputDirName,outDir,strlen(outDir));

	/*Si le fichier est une instance de la TSPLIB*/
	if(isTSPFile == true)
	{
		simpleEdge *e_list;

		e_list = readTsp((char *)filename,type);

		if( !InitGraph_from_list(e_list,tsp_n_nodes,tsp_m_edges,&G) )
		{
			cout << "ERROR:kHPPMASTER().Erreur pendant la lecture du graphe." << endl;

//			exit(Fatal);
		}

		k = k_val;

		delete []e_list;
	}
	else
	{
		/*Initialisation du graphe � partir du fichier*/
		//strncpy(File_Name,(char *)kEC_File,strlen((char *)kEC_File));

		if(!InitGraph_from_file((char *)filename,&G,&k))
		{
			cout << "ERROR:kHPPMASTER().Erreur pendant la lecture du graphe." << endl;

//			exit(Fatal);
		}
	}

    //Initialisation du graphe de demandes
	InitGraph_from_vect(&vectDemande,G.n_Nodes,vectDemande.size(),&G_dem);

	isTSPGraph = isTSPFile;
	graphType = type;

//	best_sol_num_noeud = 0;

	//*****************************************************
	//Construction des graphes orientes pour chaque demande

//	PrintGraph(&G);

	grapheDemande = new I_digraph * [vectDemande.size()];

	int s,t;

	for(int i=0;i<vectDemande.size();i++)
	{
		//Construction de la liste des arcs du graphe
		s = vectDemande[i].source;
		t = vectDemande[i].dest;

		//On les marque dans le graphe G
		G.Nodes[s-1].terminal = true;
		G.Nodes[t-1].terminal = true;

		//Construction du graphe oriente correspondant
		grapheDemande[i] = new I_digraph(s,t,&G);
//		grapheDemande[i]->printTable(cout);
	}
}



void extra_cuts::update(CUTINFO &Nc_info, IloNumArray &solX, IloArray<IloNumArray> &solF)
{
    I_Graph *g = &G;
    c_info = &Nc_info;

    for(int i=0;i < g->m_Edges;i++)
    {
		//g->Edges[i].X = solX[i];
		//g->Edges[i].back->X = solX[i];

        int node1 = g->Edges[i].back->adjac->id - 1;
        int node2 = g->Edges[i].adjac->id - 1;
//        cout << "edge " << node1  << " -- " << node2 << endl;
//        cout << "Trying to get the edge : " << c_info->G.idFromEdge( c_info->G.findEdge( c_info->G.nodeFromId( node1 ) , c_info->G.nodeFromId( node2 ) ) ) << endl;

	g->Edges[i].X = solX[ c_info->G.idFromEdge( c_info->G.findEdge( c_info->G.nodeFromId( node1 ) , c_info->G.nodeFromId( node2 ) ) ) ];
	g->Edges[i].back->X = solX[ c_info->G.idFromEdge( c_info->G.findEdge( c_info->G.nodeFromId( node1 ) , c_info->G.nodeFromId( node2 ) ) ) ];

		//cout << "Edge " << i << " " << xVal_[i] << " " << g->Edges[i].X << endl;
	}

	I_digraph *dg;
	int aux = 0;

	for(int i=0;i < vectDemande.size();i++)
	{
		dg = grapheDemande[i];    //getMaster()->getDemandGraph(i+1);

		for(int j=0;j<dg->getNbArc();j++)
		{
            		if(dg->larc[j].pedge != NIL_E)
            		{
                		dg->larc[j].cap = solX[dg->larc[j].pedge->num-1];
                		dg->larc[j].f = solX[dg->larc[j].pedge->num-1];   //solX(g->m_Edges+aux+j);
            		}
            		else
                    {
                		dg->larc[j].cap = INFINIT;
                		dg->larc[j].f = INFINIT;
            		}

			/*cout << "A[" << i << "," << j << "," <<	g->m_Edges+aux+j << "] = ";
			cout << dg->larc[j];

			if(dg->larc[j].pedge != NIL_E)
				cout << " e = " << dg->larc[j].pedge->X << " " << dg->larc[j].pedge->num;

			cout << endl;

			if((dg->larc[j].pedge != NIL_E) && (dg->larc[j].f > dg->larc[j].pedge->X))
				cout << "ERREUR" << endl;*/

			/*if(dg->larc[j].pedge != NIL_E && dg->larc[j].pedge->X > 0)
			{
				dg->larc[j].cap = dg->larc[j].pedge->X;
			}*/
		}

		aux += dg->getNbArc();
	}
}



/*****************************************************/
/*Procedure de s�paration globale pour les partition */
/*Retourne le nombre total de contraintes detectee   */
/*Type == 1 pour le normal*/
/*Type == 2 pour s�paration annexe*/
int extra_cuts::separation_globale(deque < int * > &seperated_cons, int type, FILE *sortie)
{
	//return 0;

	/*S�paration des SP-partition et Partition*/
	b_Edge *part_list = NIL_BE;
	long rhs_part;
	Bool part_result;

	I_Graph *gr,*gr_dem,gr_frac;
	I_digraph *dg;
	Bool *marquage;
	long u,v,n_e_frac,cycle_sz;
	simpleEdge *frac_edge_list;
	b_Node *frac_cycle,*bncour,*bncour2;
	I_Node *nptr;

	bool p_ok=true;
	Bool cycle_type;

	int res_p,res_sp,res_dc,res_ca,res_tpc;
	int n_p,n_sp,n_fp;

	res_p = res_sp = res_dc = res_tpc = res_ca = 0;
	n_p = 0;

	gr = &G;
	gr_dem = &G_dem;
	n = gr->n_Nodes;
	m = gr->m_Edges;

	//Separation des contraintes de sp-partition
#ifdef printings
	cout << "Separation" << endl;
#endif // printings
	//create_branch_file();

   	/*res_ca = 0;
	res_ca = AjouteCoupeAggregee(gr,k);
	getMaster()->nbCoupeAggregee += res_ca;

    	if(res_ca > 0)
        return res_ca;*/

	res_dc = 0;
	res_dc = AjouteDoubleCut(seperated_cons, gr,k,sortie);
	nbDoubleCut += res_dc;

	if(res_dc > 0)
	{
		cout << "Double cut violee - " << res_dc << " by process" << c_info->world_rank << endl;
		return res_dc;
	}

	res_tpc = 0;
	res_tpc = AjouteTriplePathCut(seperated_cons, gr,k,sortie);
	nbTriplePathCut += res_tpc;

	if(res_tpc > 0)
	{
		cout << "Triple path cut violee - " << res_tpc << " by process" << c_info->world_rank << endl;
		return res_tpc;
	}

	//create_branch_file();

	if(is_pair(k) == false) /// instead of the defined False
	{
#ifdef printings
        cout << "Contraction" << endl;
#endif // printings

		operationTheta(gr,k);

		//create_branch_file();

		//Separation des contraintes de partition
		//Le graphe reduit est constitue uniquement de terminaux.
		//On commence par prendre tous les sommets terminaux.
		//S'ils sont en nombre pair, on contracte
		//deux terminaux qui ne sont pas dans la meme demande
		//et qui minimise le poids total de la nouvelle partition.
		frac_cycle = NIL_BN;
		nptr = gr->Nodes[0].n_sh;
		cycle_sz = 0;

		do
		{
			if(nptr->terminal)
			{
				Insere_Node(&frac_cycle,nptr->id);
				cycle_sz++;

				//cout << "U = " << nptr->id << endl;
			}

			nptr = nptr->next_gr;

		}while(nptr != gr->Nodes[0].n_sh);

		/*cout << "FRAC CYCLE 1 = " << cycle_sz << endl;
		Print_b_Node_Set(frac_cycle,stdout);
		cout << endl;*/

		if(is_pair(cycle_sz) == true) /// instead of the defined True
		{
			//cout << "Cycle Pair = " << cycle_sz << endl;

			//On cherche deux sommets consecutifs u et v tel que
			//u et v ne forment pas une demande et (k-1)/2 < x([u,v]) <= k/2.
			bool contraction = false;
			bncour = frac_cycle;
			while((bncour != NIL_BN) && (!contraction))
			{
				bncour2 = frac_cycle;
				while((bncour2 != NIL_BN) && (!contraction))
				{
					if(bncour2->id != bncour->id)
					{
						//On essaie de contracter deux terminaux
						//qui ne forment pas la meme demande
						/*int i=0;
						while((i<vectDemande.size()) &&
						!(
						((getMaster()->getDemande(i).source == bncour->id) &&
						(getMaster()->getDemande(i).dest == bncour2->id) )
						||
						((getMaster()->getDemande(i).source == bncour2->id) &&
						(getMaster()->getDemande(i).dest == bncour->id) )
						))
						i++;*/


						/************************************/
						/*  On verifie si les deux sommets  */
						/************************************/

						b_Node *w_dem = NIL_BN;

						Insere_Node(&w_dem,bncour->id);
						Insere_Node(&w_dem,bncour2->id);

						/*cout << "A contracter ?" << endl;
						Print_b_Node_Set(w_dem,stdout);*/

						if(admissible(gr_dem,w_dem))
						{
							b_Edge *edge_u_v,*becour;
							edge_u_v = get_edge_u_v(gr,bncour->id,bncour2->id,NO_0_EDGES,CAP_VALUE_X);

							//Print_b_Edge_Set(edge_u_v,stdout);

							becour = edge_u_v;
							double sommeX = 0;
							while(becour != NIL_BE)
							{
								sommeX += becour->cap;

								becour = becour->next;
							}

							//cout << "SOMME X = " << sommeX << endl;

							Delete_b_Edge_Set(&edge_u_v);

							if(((k-1.0)/2 < sommeX) && (sommeX <= (1.0*k)/2))
							{
								save_contraction_info(gr,1);

								//On contracte bncour et bncour2
								b_Node *W = NIL_BN;

								Insere_Node(&W,bncour->id);
								Insere_Node(&W,bncour2->id);

								contract_set_w(gr,W);

								//Print_b_Node_Set(W,stdout);

								Delete_b_Node_Set(&W);

								contraction = true;
							}

						}

						Delete_b_Node_Set(&w_dem);
					}

					bncour2 = bncour2->next;
				}

				bncour = bncour->next;
			}

			Delete_b_Node_Set(&frac_cycle);

			res_p = 0;

			if(contraction)
			{
				//S'il y a eu contraction on reconstruit
				//la liste des sommets

				frac_cycle = NIL_BN;
				nptr = gr->Nodes[0].n_sh;
				cycle_sz = 0;

				do
				{
					if(nptr->terminal)
					{
						Insere_Node(&frac_cycle,nptr->id);
						cycle_sz++;

						//cout << "U = " << nptr->id << endl;
					}

					nptr = nptr->next_gr;

				}while(nptr != gr->Nodes[0].n_sh);

				/*cout << "FRAC CYCLE 2 = " << cycle_sz << endl;
				Print_b_Node_Set(frac_cycle,stdout);
				cout << endl;*/

				//create_branch_file();

				res_p = 0;
				res_p = AjoutePartition(seperated_cons, gr,k,frac_cycle,cycle_sz,sortie);
				nbPartition += res_p;

				recall_contraction_info(gr,1);
			}
		}
		else
		{
			res_p = 0;
			res_p = AjoutePartition(seperated_cons, gr,k,frac_cycle,cycle_sz,sortie);
			nbPartition += res_p;
		}

		Delete_b_Node_Set(&frac_cycle);

		if(res_p>0)
		{
			cout << "Partition violee - " << res_p << " by process" << c_info->world_rank << endl;
			return res_p;
		}

		//create_branch_file();


		res_sp = 0;
		res_sp = AjouteSPPartition2(seperated_cons, gr,k,frac_cycle,cycle_sz,sortie);

		//res_sp = AjouteSPPartition(gr,k,frac_cycle,cycle_sz,sortie);
		//res_sp = AjouteSPPartition2(gr,k,frac_cycle,cycle_sz,sortie);

		nbSPPartition += res_sp;

		if(res_sp>0)
		{
			cout << "SP-Partition violee - " << res_sp << " by process " << c_info->world_rank << endl;
			return res_sp;
		}
	}
	return (res_p + res_sp + res_dc + res_tpc);
}


int extra_cuts::AjouteFPartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie)
{
/*	long *f_part_list=NULL;
	long *V0_F,V0_F_sz;

	long rhs_f_part;
	Bool f_part_result;
	int nb_f_part = 0;
	bool fp_ok;

	fp_ok = true;

	/*S�paration des F-Partition*/
/*	f_part_result = False;
	f_part_result = separation_f_partition(gr,k,frac_cycle,cycle_sz,&f_part_list,&V0_F,&V0_F_sz,&rhs_f_part,sortie);


	if(f_part_result == True)
	{
		/*Ajout de la contrainte dans le pool*/
/*		ABA_BUFFER<ABA_CONSTRAINT *> fp_constraints(master_,1);
		F_PARTITION *new_fp;

		new_fp = new F_PARTITION(master_,f_part_list,cycle_sz+1,V0_F,V0_F_sz,rhs_f_part);
		fp_constraints.push(new_fp);

		int n_fp = addCons(fp_constraints,get_master()->get_pool_f_partition(),0,0);
		//addCons(fp_constraints);

		nb_f_part++;


		if(n_fp != 1)
		{
			fp_ok = false;
		}
		else
		{
			fp_ok = true;
		}

		//cout << endl << "F_Partition viol�e" << endl;
		//Print_b_Edge_Set(f_part_list,stdout);
		//cout << "Orginal" << endl;
		//Print_b_Edge_Set_Gr(gr,f_part_list,stdout);
		//create_branch_file(Te,m);

/*	}
	else
	{
		nb_f_part += 0;
		//cout << endl << "Pas de F-Partition viol�e" << endl;
	}


	if(f_part_result == True)
	{
		if(fp_ok)
		{
			get_master()->add_nb_f_part(nb_f_part);
			return nb_f_part;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		return 0;
	}*/

	return 0;
}



int extra_cuts::AjoutePartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie)
{
	long *part_list=NULL;
	long rhs_part;
	Bool part_result;
	int nb_part = 0;
	bool p_ok;

	p_ok = true;

	part_result = false; /// instead of the defined False
	part_result = separation_partition_steiner(gr,k,frac_cycle,cycle_sz,&part_list,&rhs_part,sortie);

	if(part_result == true) /// instead of the defined True
	{
		/*Ajout de la contrainte dans le pool*/
//		ABA_BUFFER<ABA_CONSTRAINT *> p_constraints(master_,1);
//		PARTITION *new_p;
//
//		new_p = new PARTITION(master_,part_list,cycle_sz+1,rhs_part);
//		p_constraints.push(new_p);
//
//		int n_p = addCons(p_constraints,getMaster()->getPoolPartition(),0,0);

        vector <int> *_cons;
        double aCoeff = 0;

        _cons = new vector <int> ();

        _cons->push_back(0);
        _cons->push_back(P);
        (*_cons)[0] += 2;

        PARTITION *new_p;
        new_p = new PARTITION(part_list,cycle_sz+1,rhs_part);

        for(int i=0; i < gr->m_Edges; i++)
        {

            aCoeff = (int) new_p->coeff(*c_info, i);

            if (aCoeff > 0)
            {
                _cons->push_back(aCoeff);
                _cons->push_back(i);
                (*_cons)[0] += 2;
            }
            //aCtrExpr += aCoeff * c_info->xSol[i];
        }

        _cons->push_back(GET);
        _cons->push_back(rhs_part);
        (*_cons)[0] += 2;

        seperated_cons.push_back( & ( (*_cons)[0] ) );
        int n_p = 1;

		//addCons(p_constraints);


		nb_part++;

		if(n_p != 1)
		{
			p_ok = false;
		}
		else
		{
			p_ok = true;
		}

		//cout << endl << "Partition viol�e" << endl;
		//Print_b_Edge_Set(part_list,sortie_frac_gk);

		//cout << endl << "Partition viol�e" << endl;
		//Print_b_Edge_Set(part_list,stdout);
		//cout << "Orginal" << endl;
		//Print_b_Edge_Set_Gr(gr,part_list,stdout);
		//create_branch_file(Te,m);
	}
	else
	{
		nb_part += 0;
		//cout << endl << "Pas de Partition viol�e" << endl;
	}

	if(part_result == true)   /// instead of defined True
	{
		if(p_ok)
		{
			//getMaster()->add_nb_part(nb_part);
			return nb_part;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		return 0;
	}

	return 0;
}



int extra_cuts::AjouteDoubleCut(deque < int * > &seperated_cons, I_Graph *gr,int k,FILE *sortie)
{
/// To verify if this can be done inside the next (parallel) loop
	I_digraph *dg;

	for(int i=0;i < vectDemande.size();i++)
	{

		dg = grapheDemande[i];

		//cout << "DEM == " << getMaster()->getDemande(i+1).source << " " << getMaster()->getDemande(i+1).dest
		//<< endl;

		for(int j=0;j < dg->getNbArc();j++)
		{
            if(dg->larc[j].pedge != NIL_E && dg->larc[j].pedge->X > 0)
			{
				dg->larc[j].cap = dg->larc[j].pedge->X;

				//cout << "A = " << dg->larc[j] << endl;
			}
		}
	}

	//create_branch_file();

	int n,m;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	int nbDC = 0;

	//Pour chaque demande, on essaie de trouver une contrainte
//	#pragma omp parallel for private(dg)
	for(int i=0;i < vectDemande.size();i++)
	{
	    I_Graph *_gr;
        _gr = new I_Graph;
	    copy_Graph(_gr,gr);

		int s = vectDemande[i].source;
		int t = vectDemande[i].dest;
		int t2;

		/*I_digraph **/
		dg = grapheDemande[i];

		b_Node *pi[5];
		b_Edge *delta_s,*delta_t2,*ecour,*ecour2;

		//On recherche un autre terminal t1 candidat

		/*Insere_Node(&pi[0],s);
		Insere_Node(&pi[4],t);*/

		/*cout << "s1 == " << s << endl;
		cout << "t1 == " << t << endl;*/

		delta_s = get_delta_v(_gr,s,NO_0_EDGES,CAP_VALUE_X);
		double X_D;

		X_D = 0;
		ecour = delta_s;
		while(ecour != NIL_BE)
		{
			X_D += ecour->cap;

			ecour = ecour->next;
		}

		//cout << "XD == " << X_D << endl;

		if(X_D == (double)k)
		{

			//On cherche un voisin de s qui induit une coupe serree
			ecour = delta_s;
			while(ecour != NIL_BE)
			{
				bool candidat_ok = false;
				int candidat;

				candidat = ecour->node2;

				if((candidat != t) && (_gr->Nodes[candidat-1].terminal == true))
				{
					delta_t2 = get_delta_v(_gr,candidat,NO_0_EDGES,CAP_VALUE_X);
					ecour2 = delta_t2;

					X_D = 0;
					while(ecour2 != NIL_BE)
					{
						X_D += ecour2->cap;

						ecour2 = ecour2->next;
					}

					Delete_b_Edge_Set(&delta_t2);

					if(X_D == (double)k)
						candidat_ok = true;
				}

				if(candidat_ok)
				{
					//t2 = ecour->node2;
					t2 = candidat;

					//Insere_Node(&pi[1],t2);

					//cout << "t2 == " << t2 << endl;

					//Maintenant, on cherche une L-st-path-cut
					//avec V1 = {t2}. Pour ca, on utilise le graphe
					//oriente.
					//Il faut d'abord modifier le graphe oriente
					//pour preciser qu'on cherche une coupe contenant
					//les aretes de delta(s)\st2.
					//On met 0 comme capacite pour les arcs de
					//delta+(s)\(s,t2) et capacite INF pour l'arc ()

					I_diarc *a;
					list<int> num_arc_modif;
					list<double> cap_arc_modif;

					a = dg->source.premier;
					while(a != NIL_A)
					{
						if(a->cap > 0)
						{
							//cout << "Modif A == " << *a << endl;

							cap_arc_modif.push_back(a->cap);
							num_arc_modif.push_back(a->num);

							if(a->dest->num_reel != t2)
								a->cap = 0;
							else
								a->cap = INFINI;
						}

						a = a->suivant;
					}

					//create_branch_file();

					//Ensuite on calcule un flot max entre s et t
					list<I_dinode *> W;
					double flotMax = dg->calculFlot(W);

					/*cout << "Flot == " << flotMax << endl;
					cout << "Ensemble W == " << endl;

					list<I_dinode *>::iterator it;
					for(it=W.begin();it!=W.end();it++)
					{
						cout << *(*it) << " ";
					}
					cout << endl;*/


					//On reinitialise les arcs modifies
					list<int>::iterator itl_num;
					list<double>::iterator itl_cap = cap_arc_modif.begin();

					for(itl_num = num_arc_modif.begin();itl_num != num_arc_modif.end();itl_num++)
					{
						dg->larc[*itl_num-1].cap = *itl_cap;
						//cout << "Restore = " << dg->larc[*itl_num-1] << endl;

						itl_cap++;
					}

					//Constitution de la partition.
					//Pour cela, on calcule les distances
					//entre s et tous les sommets
					//Avant cela, il faut supprimer les aretes
					//correspondant a la coupe delta+(W)
					list<I_diarc *> delta_w;
					list<I_diarc *>::iterator ita;

					dg->getDelta(W,delta_w);

					b_Edge *edge_modif,*ecour3;

					edge_modif = NIL_BE;

					for(ita = delta_w.begin();ita != delta_w.end();ita++)
					{
						//cout << "A = " << *(*ita) << endl;

						if((*ita)->cap > 0)
						{
							//cout << "A = " << *(*ita) << endl;

							if((*ita)->pedge != NIL_E)
							{
								ecour3 = new b_Edge;
								ecour3->node1 = (*ita)->pedge->back->adjac->id;
								ecour3->node2 = (*ita)->pedge->adjac->id;
								ecour3->cap = (*ita)->pedge->X;
								ecour3->num = (*ita)->pedge->num;
								ecour3->next = NIL_BE;

								Insere_Edge(&edge_modif,ecour3);

//								(*ita)->pedge->X = 0;
//								(*ita)->pedge->back->X = 0;

								_gr->Edges[(*ita)->pedge->num-1].X = 0;
								_gr->Edges[(*ita)->pedge->num-1].back->X = 0;
							}
							else
							{
								cout << "Erreur dans la coupe" << endl;

								//Dans ce cas, s'il s'agit d'un
								//arc de la forme (u',u""), alors
								//l'arc (s,u') a une capacite +INF
								//Il s'en suit que l'arc (u',t) peut
								//etre mis dans la coupe sans changer le poids
								//On marque donc le pedge correspondant

								/*I_dinode *u_pri,*u_sec;
								I_diarc *a_sec;

								u_pri = (*ita)->source;
								u_sec = (*ita)->dest;

								a_sec = u_sec->premier;
								while(a_sec != NIL_A && a_sec->dest->num_reel != t+3*n)
									a_sec = a_sec->suivant;

								ecour3 = new b_Edge;
								ecour3->node1 = a_sec->pedge->back->adjac->id;
								ecour3->node2 = a_sec->pedge->adjac->id;
								ecour3->cap = a_sec->pedge->X;
								ecour3->num = a_sec->pedge->num;
								ecour3->next = NIL_BE;

								Insere_Edge(&edge_modif,ecour3);

//								a_sec->pedge->X = 0;
//								a_sec->pedge->back->X = 0;
                                _gr->Edges[a_sec->pedge->num-1].X = 0;
								_gr->Edges[a_sec->pedge->num-1].back->X = 0;*/
							}
						}
					}

					//Print_b_Edge_Set(edge_modif,stdout);

					double *T = NULL;
					int *P = NULL;
					int *PredArete = NULL;
					double *L = NULL;
					Dijkstra(_gr,s,L,P,PredArete,CAP_VALUE_1);

					/*cout << "Dijkstra" << endl;
					for(int j=0;j<n;j++)
					{
						cout << "L[" << j << "] = " << L[j] << endl;
					}

					for(int j=0;j<n;j++)
					{
						cout << "P[" << j << "] = " << P[j] << endl;
					}*/


					//Il ne faut pas oublier de remettre
					//les valeurs des aretes modifiees

					ecour3 = edge_modif;
					while(ecour3 != NIL_BE)
					{
						_gr->Edges[ecour3->num-1].X = ecour3->cap;
						_gr->Edges[_gr->m_Edges+ecour3->num-1].X = ecour3->cap;

						ecour3 = ecour3->next;
					}

					Delete_b_Edge_Set(&edge_modif);

					//Il faut verifier qu'il s'agit d'une L-path-cut
					//Pour cela, il faut que L[t-1] < INFINI

					if(L[t-1] < INFINI)
					{
						for(int j=0;j<5;j++)
							pi[j] = NIL_BN;

						//Formule: iVAlter =  iVAlter + par
						//par = -1*par
						int iVAlternatif;	//Vaudra 1 ou 3
						int par = 1;		//Vaudra -1 ou 1

						for(int j=0;j<n;j++)
						{
							/*int l;

							if(((int)L[j]) <= 3)
								l = (int)L[j];
							else
								l = 4;*/

							int ld;

							if(j == s-1)
								ld = 0;
							else if(j == t-1)
								ld = 4;
							else
							{
								if(((int)L[j]) <= 3)
									ld = (int)L[j];
								else
								{
									//Un sommet autre que t1 qui est à distance >= 4 est mis dans V3
									//Ceux qui sont à distance infinie seront mis dans V1 et V3:
									//une partie dans V1 et l'autre dans V3
									if(((int)L[j]) < INFINI)
										ld = 3;
									else
									{
										par = -1*par;
										iVAlternatif = 2 + par;

										ld = iVAlternatif;
									}
								}
							}

							/*cout << "J == " << j << " L = " << L[j] << endl;
							cout << "l == " << l << endl;*/

							Insere_Node(&pi[ld],j+1);
						}

						//Delete_b_Edge_Set(&edge_modif);

						int *partition = NULL;
						vector<int> F;
						int rhs;

						Bool dc_result = false; /// instead of the defined False
						dc_result = separation_double_cut(_gr,pi,k,s,t,partition,F,rhs,sortie);

						//static int aDc = 1;

						if(dc_result == true /*&& aDc <= 200*/) ///instead of the defined True
						{
                            //aDc++;

                            /////////////////////////////////////////////////////////////////////

                            vector <int> *_cons;
                            double aCoeff = 0;

                            _cons = new vector <int> ();

                            _cons->push_back(0);
                            _cons->push_back(DC);
                            (*_cons)[0] += 2;

                            DOUBLE_CUT *new_dc;
                            new_dc = new DOUBLE_CUT(partition,F,rhs);

                            for(int i=0; i < gr->m_Edges; i++)
                            {

                                aCoeff = new_dc->coeff(*c_info, i);

                                if (aCoeff > 0)
                                {
                                    //cout << "Here X" << i << " coef = " << aCoeff << endl;
                                    _cons->push_back(aCoeff);
                                    _cons->push_back(i);
                                    (*_cons)[0] += 2;
                                }
                                //aCtrExpr += aCoeff * c_info->xSol[i];
                            }

                            _cons->push_back(GET);
                            _cons->push_back(rhs);
                            (*_cons)[0] += 2;

                            seperated_cons.push_back( & ( (*_cons)[0] ) );

                            nbDC++;

                            /////////////////////////////////////////////////////////////////////


							//Ajout de la contrainte dans le pool
//							ABA_BUFFER<ABA_CONSTRAINT *> p_constraints(master_,1);
//							DOUBLE_CUT *new_dc;
//
//							new_dc = new DOUBLE_CUT(getMaster(),partition,F,rhs);
//							p_constraints.push(new_dc);
//
//                            l1.lock();
//							addCons(p_constraints,getMaster()->getPoolDoubleCut(),0,0);
//							//addCons(p_constraints);
//
//							//cout << "DC violee" << endl << endl;
//
//							nbDC++;
//							l1.unlock();
						}
						else
						{
							//cout << endl << "Pas de Partition viol�e" << endl;
						}

						//On supprime l'ensemble pi qui a ete cree
						for(int j=0;j<5;j++)
							Delete_b_Node_Set(&pi[j]);
					}
					else
					{
#ifdef printings
						cout << "T injoignable" << endl;
#endif // printings
					}

					delete[] L;
					delete[] P;
					delete[] PredArete;
					delete[] T;

				}
				else
				{
					//cout << "Mauvais candidat = " << ecour->node2 << endl;
				}

				ecour = ecour->next;
			}
		}

		Delete_b_Edge_Set(&delta_s);
		DeleteGraph(_gr);
	}

	return nbDC;
}



int extra_cuts::AjouteTriplePathCut(deque < int * > &seperated_cons, I_Graph *gr,int k,FILE *sortie)
{
	I_digraph *dg;

	for(int i=0;i < vectDemande.size();i++)
	{
		dg = grapheDemande[i];

		for(int j=0;j<dg->getNbArc();j++)
		{
			if(dg->larc[j].pedge != NIL_E && dg->larc[j].pedge->X > 0)
			{
				dg->larc[j].cap = dg->larc[j].pedge->X;
			}
		}
	}
	//create_branch_file();

	int *partition = NULL;
	int n,m;


	n = gr->n_Nodes;
	m = gr->m_Edges;

	int nbTPC = 0;

	I_Graph *gr_dem = &G_dem;

	//Pour chaque source s, on regarde si s est associée à au moins 2 destinations  t1,t2
	//Si c'est le cas, on essaie de construire une triple path cut
	//Enfin, on vérifie si la contrainte est violée ou non

//    #pragma parallel for private(dg)
	for(int i=0;i<source.size();i++)
	{
        I_Graph *_gr;
        _gr = new I_Graph;
	    copy_Graph(_gr,gr);

	    int s,t1,t2;

		s = source[i];

		//On calcule les destinations qui sont associées à t1, t2
		b_Edge *delta_w_dem = get_delta_v(gr_dem,s,NO_0_EDGES,CAP_VALUE_1);

		vector<pair<int,int> > listDest;	//Le premier élément de la pair est le sommet
											//le deuxième élément est le numéro de la demande
											//qui est donné par le numéro de l'arete dans le graphe de demande
		b_Edge *be_cour;

		be_cour = delta_w_dem;
		while(be_cour != NIL_BE)
		{
			listDest.push_back(make_pair<int,int>(be_cour->node2,be_cour->num));

			be_cour = be_cour->next;
		}

		Delete_b_Edge_Set(&delta_w_dem);

		if(listDest.size() >= 2)
		{
			/*cout << "Source candidate = " << s << endl;
			for(int j=0;j<listDest.size();j++)
			{
				cout << " T = " << listDest[j].first << " " << listDest[j].second << endl;
			}*/

			b_Node *pi[6];

			for(int j=0;j<listDest.size()-1;j++)
			{
				t1 = listDest[j].first;

				for(int r=j+1;r<listDest.size();r++)
				{
					t2 = listDest[r].first;

					//cout << "S = " << s << " T1 = " << t1 << " T2 = " << t2 << endl;

					dg = grapheDemande[listDest[j].second - 1];
					//cout << dg << endl;

					//On calcule maintenant la 2e L-path-cut (V0,V1 U V4,V2,V3,V5)
					//Grace à elle, on peut déterminer toute la partition


					I_diarc *a;
					list<int> num_arc_modif;
					list<double> cap_arc_modif;

					a = dg->source.premier;
					while(a != NIL_A)
					{
						if(a->cap > 0)
						{
							//cout << "Modif A == " << *a << endl;

							if(a->dest->num_reel != t2)
							{
								/*a->cap = 0;
								cap_arc_modif.push_back(a->cap);
								num_arc_modif.push_back(a->num);*/
							}
							else
							{
								cap_arc_modif.push_back(a->cap);
								num_arc_modif.push_back(a->num);

								a->cap = INFINI;
							}
						}

						a = a->suivant;
					}

					//create_branch_file();

					//Ensuite on calcule un flot max entre s1 et t1
					list<I_dinode *> W;
					double flotMax = dg->calculFlot(W);

					/*cout << "Flot == " << flotMax << endl;
					cout << "Ensemble W == " << endl;

					list<I_dinode *>::iterator it;
					for(it=W.begin();it!=W.end();it++)
					{
						cout << *(*it) << " ";
					}
					cout << endl;*/


					//On reinitialise les arcs modifies
					list<int>::iterator itl_num;
					list<double>::iterator itl_cap = cap_arc_modif.begin();

					for(itl_num = num_arc_modif.begin();itl_num != num_arc_modif.end();itl_num++)
					{
						dg->larc[*itl_num-1].cap = *itl_cap;
						//cout << "Restore = " << dg->larc[*itl_num-1] << endl;

						itl_cap++;
					}

					//Constitution de la partition.
					//Pour cela, on calcule les distances
					//entre s et tous les sommets
					//Avant cela, il faut supprimer les aretes
					//correspondant a la coupe delta+(W)
					list<I_diarc *> delta_w;
					list<I_diarc *>::iterator ita;
					double poidsLPathCut = 0;		//Poids de la L-path-cut

					dg->getDelta(W,delta_w);

					b_Edge *edge_modif,*ecour3;

					edge_modif = NIL_BE;

					for(ita = delta_w.begin();ita != delta_w.end();ita++)
					{
						//cout << "A = " << *(*ita) << endl;

						if((*ita)->cap > 0)
						{
							//cout << "A = " << *(*ita) << endl;

							if((*ita)->pedge != NIL_E)
							{
								ecour3 = new b_Edge;
								ecour3->node1 = (*ita)->pedge->back->adjac->id;
								ecour3->node2 = (*ita)->pedge->adjac->id;
								ecour3->cap = (*ita)->pedge->X;
								ecour3->num = (*ita)->pedge->num;
								ecour3->next = NIL_BE;

								Insere_Edge(&edge_modif,ecour3);

								poidsLPathCut += (*ita)->pedge->X;

//								(*ita)->pedge->X = 0;
//								(*ita)->pedge->back->X = 0;

								_gr->Edges[(*ita)->pedge->num-1].X = 0;
								_gr->Edges[(*ita)->pedge->num-1].back->X = 0;
							}
							else
							{
								cout << "Erreur dans la coupe" << endl;

								//Dans ce cas, s'il s'agit d'un
								//arc de la forme (u',u""), alors
								//l'arc (s,u') a une capacite +INF
								//Il s'en suit que l'arc (u',t) peut
								//etre mis dans la coupe sans changer le poids
								//On marque donc le pedge correspondant

								/*I_dinode *u_pri,*u_sec;
								I_diarc *a_sec;

								u_pri = (*ita)->source;
								u_sec = (*ita)->dest;

								a_sec = u_sec->premier;
								while(a_sec != NIL_A && a_sec->dest->num_reel != t+3*n)
									a_sec = a_sec->suivant;

								ecour3 = new b_Edge;
								ecour3->node1 = a_sec->pedge->back->adjac->id;
								ecour3->node2 = a_sec->pedge->adjac->id;
								ecour3->cap = a_sec->pedge->X;
								ecour3->num = a_sec->pedge->num;
								ecour3->next = NIL_BE;

								Insere_Edge(&edge_modif,ecour3);

								a_sec->pedge->X = 0;
								a_sec->pedge->back->X = 0;*/
							}
						}
					}

					//Print_b_Edge_Set(edge_modif,stdout);

					double *T = NULL;
					int *P = NULL;
					int *PredArete = NULL;
					double *L = NULL;
					Dijkstra(_gr,s,L,P,PredArete,CAP_VALUE_1);

					//Il ne faut pas oublier de remettre
					//les valeurs des aretes modifiees

					ecour3 = edge_modif;
					while(ecour3 != NIL_BE)
					{
						_gr->Edges[ecour3->num-1].X = ecour3->cap;
						_gr->Edges[_gr->m_Edges+ecour3->num-1].X = ecour3->cap;

						ecour3 = ecour3->next;
					}

					Delete_b_Edge_Set(&edge_modif);

					//Il faut verifier qu'il s'agit d'une L-path-cut
					//Pour cela, il faut que L[t1-1] < INFINI

					if(L[t1-1] < INFINI)
					{
						//cout << "L-path-cut trouvee - Poids = " << poidsLPathCut << endl;

						for(int j=0;j<6;j++)
							pi[j] = NIL_BN;

						int iV1V3 = 1;	//Vaudra 1 ou 3
						int par = 1;	//Vaudra 0 ou 1

						for(int j=0;j<n;j++)
						{
							int ld;

							if(j == s-1)
								ld = 0;
							else if(j == t1-1)
								ld = 5;
							else if(j == t2-1)
								ld = 4;
							else
							{
								if(((int)L[j]) <= 3)
									ld = (int)L[j];
								else
								{
									//Un sommet autre que t1 qui est à distance >= 4 est mis dans V3
									//Ceux qui sont à distance infinie seront mis dans V1 et V3:
									//une partie dans V1 et l'autre dans V3
									if(((int)L[j]) < INFINI)
										ld = 3;
									else
									{
										par = (par+1)%2;
										iV1V3 = 1 + 2*par;

										ld = iV1V3;
									}
								}
							}

							/*cout << "J == " << j << " L = " << L[j] << endl;
							cout << "l == " << l << endl;*/

							Insere_Node(&pi[ld],j+1);
						}

						/*for(int j=0;j<6;j++)
						{
							cout << "V" << j << " = ";
							Print_b_Node_Set(pi[j],stdout);
							cout << endl;
						}*/


						int *partition = NULL;
						vector<int> F;
						int rhs;

						Bool tpc_result = false;  /// instead of defined False
						tpc_result = separation_triple_path_cut(_gr,pi,k,s,t1,t2,partition,F,rhs,sortie);

						static int aTpc = 1;

						if(tpc_result == true && aTpc <= 200) /// instead of defined True
						{
							aTpc++;

							/////////////////////////////////////////////////////////////////////

                            vector <int> *_cons;
                            double aCoeff = 0;

                            _cons = new vector <int> ();

                            _cons->push_back(0);
                            _cons->push_back(TPC);
                            (*_cons)[0] += 2;

                            TRIPLE_PATH_CUT *new_tpc;
                            new_tpc = new TRIPLE_PATH_CUT(partition,F,rhs);

                            for(int i=0; i < gr->m_Edges; i++)
                            {

                                aCoeff = (int) new_tpc->coeff(*c_info, i);

                                if (aCoeff > 0)
                                {
                                    _cons->push_back(aCoeff);
                                    _cons->push_back(i);
                                    (*_cons)[0] += 2;
                                }
                                //aCtrExpr += aCoeff * c_info->xSol[i];
                            }

                            _cons->push_back(GET);
                            _cons->push_back(rhs);
                            (*_cons)[0] += 2;

                            seperated_cons.push_back( & ( (*_cons)[0] ) );

                            nbTPC++;

                            /////////////////////////////////////////////////////////////////////


							//Ajout de la contrainte dans le pool
//							ABA_BUFFER<ABA_CONSTRAINT *> p_constraints(master_,1);
//							TRIPLE_PATH_CUT *new_tpc;
//
//							new_tpc = new TRIPLE_PATH_CUT(getMaster(),partition,F,rhs);
//							p_constraints.push(new_tpc);
//
//                            l1.lock();
//							addCons(p_constraints,getMaster()->getPoolTriplePathCut(),0,0);
//							//addCons(p_constraints);
//
//							//cout << "TPC violee" << endl << endl;
//
//							nbTPC++;
//							l1.unlock();
						}
						else
						{
							//cout << endl << "Pas de Partition violée" << endl;
						}

						//On supprime l'ensemble pi qui a ete cree
						for(int j=0;j<6;j++)
							Delete_b_Node_Set(&pi[j]);
					}
					else
					{
						//cout << "T injoignable" << endl;
					}

					delete[] L;
					delete[] P;
					delete[] PredArete;
					delete[] T;
				}
			}
		}
        DeleteGraph(_gr);
	}


	return nbTPC;
}



int extra_cuts::AjouteCoupeAggregee(deque < int * > &seperated_cons, I_Graph *gr,int k)
{
    /***************************************************/
	/*   On fait la separation des coupes aggregees    */
	/***************************************************/

	I_digraph *dg;

	int res_ca = 0;

   	int n = gr->n_Nodes;
    int nbDem = vectDemande.size();

    //Sera utilisé comme le parametre dans les labels pour
    //distinguer les numéros des arcs des différents graphes dg
    int parDG = grapheDemande[1]->getNbArc(); //On suppose que tous les graphes dg ont le m�me nombre d'arcs

    //Nombre total d'arcs
	int nbArcDG = 0;
	for(int i=0;i<vectDemande.size();i++)
        nbArcDG += grapheDemande[i]->getNbArc();

	//Construction de la liste des arcs fractionnaires et de leurs
	//correspondance
	vector<int> labelArcFracH(nbArcDG);     //un tableau qui permet d'indicer tous les arcs de tous les graphes dg
	vector<I_diarc *> *labelNumArcFracH = new vector<I_diarc *>[nbDem];
	vector<simpleEdge> *labelNumArcFracH2 = new vector<simpleEdge>[nbDem];
	vector<I_diarc *> labelNumArcFracH3;
	int nbArcFrac,nbArcFracTotal;
	vector<simpleEdge> vectAreteH;

	vector< list<I_dinode *> > vectWH;
	vector< list<I_diarc *> > vectLEntier;
	vector<int> vectDemandeWLEntier;    //Sert � stocker le num�ro de la demande � laquel le W et le LEntier sont associ�s

    int parDGFrac = 0;
    nbArcFracTotal = 0;

    //#pragma parallel for reduction(+:parDGFrac,nbArcFracTotal)
	for(int j=0;j<nbDem;j++)
	{
        nbArcFrac = 0;
        dg = grapheDemande[j];

        for(int i=0;i<dg->getNbArc();i++)
        {
            if((dg->larc[i].cap > 0) && (dg->larc[i].cap < 1))
            {
                //cout << "J = " << j << " I = " << i << endl;

                labelNumArcFracH[j].push_back(&dg->larc[i]);

                /****************************************/
                simpleEdge se;

                se.node1 = dg->larc[i].source->num_reel;
                se.node2 = dg->larc[i].dest->num_reel;
                se.cap = dg->larc[i].cap;

                /*cout << "Se1 == " << se.node1 << endl;
                cout << "Se2 == " << se.node2 << endl;*/

                labelNumArcFracH2[j].push_back(se);
                labelNumArcFracH3.push_back(&dg->larc[i]);
                /****************************************/

                nbArcFrac++;
                labelArcFracH[j*parDG+dg->larc[i].num-1] = parDGFrac+nbArcFrac;
            }
        }

        parDGFrac += labelNumArcFracH[j].size();
        nbArcFracTotal += labelNumArcFracH[j].size();
    }

    //create_branch_file();

	//On commence par mettre dans la liste d'arete de H
	//les aretes reliant deux arcs associes a la meme arete
    parDGFrac = 0;

    //Creating draft variables to create parallelism
    //Section instead of deleting manually draftvariables
    {
        vector < vector <simpleEdge> > _vectAreteHs;
        int _parDGFracs[nbDem];
        _parDGFracs[0] = 0;
        _vectAreteHs.push_back(vector <simpleEdge> ());
        for(int r=1; r<nbDem; r++)
        {
            _vectAreteHs.push_back(vector <simpleEdge> ());
            _parDGFracs[r] = _parDGFracs[r-1] + labelNumArcFracH[r-1].size();
        }

//        #pragma omp parallel for
        for(int r=0; r<nbDem; r++)
        {
            //cout << "R = " << r << endl;

            for(int i=0; i<labelNumArcFracH[r].size(); i++)
            {
                //cout << " I = " << i << " " << *labelNumArcFracH[r][i] << endl;

                //On le compare avec tous ceux qui sont dans la m�me demande que lui
                for(int j=i+1; j<labelNumArcFracH[r].size(); j++)
                {
                    //cout << "  J = " << j << endl;

                    //On identifiera ces arcs grace a leur
                    //capacite = 0.5
                    if(corresponde(labelNumArcFracH[r][i],labelNumArcFracH[r][j],n))
                    {
                        /*cout << "ac = " << *labelNumArcFracH[i] << endl;
                        cout << "ac2 = " << *labelNumArcFracH[j] << endl;*/

                        simpleEdge se;
                        se.node1 = _parDGFracs[r]+i+1;
                        se.node2 = _parDGFracs[r]+j+1;
                        se.cap = 0.5;

                        /*cout << "Se1-1 == " << se.node1 << endl;
                        cout << "Se2-1 == " << se.node2 << endl;*/
                        _vectAreteHs[r].push_back(se);
                    }
                }

                //puis avec les arc fractionnaires des autres demandes
                int parDGFrac2 = _parDGFracs[r];

                for(int r2=r+1; r2<nbDem; r2++)
                {
                    parDGFrac2 += labelNumArcFracH[r2-1].size();

                    for(int j=0; j<labelNumArcFracH[r2].size(); j++)
                    {
                        if(corresponde(labelNumArcFracH[r][i],labelNumArcFracH[r2][j],n))
                        {
                            simpleEdge se;
                            se.node1 = _parDGFracs[r]+i+1;
                            se.node2 = parDGFrac2+j+1;
                            se.cap = 0.5;

                            _vectAreteHs[r].push_back(se);
                        }
                    }
                }
            }
        }

        for(int r=0; r<nbDem; r++)
            for (int rr=0; rr<_vectAreteHs[r].size(); rr++)
                vectAreteH.push_back(_vectAreteHs[r][rr]);
    }
    //Nombre d'aretes rouge. Ca donne l'index de depart
	//des aretes noires dans la liste des aretes de H
	int pRouge = vectAreteH.size();

    //On ajoute les arcs noires
    parDGFrac = 0;

    //Creating draft variables to create parallelism
    //Section instead of deleting manually draftvariables
    {
        vector <simpleEdge> _vectAreteHs[vectDemande.size()];
        vector < list<I_dinode *> >  _vectWHs[vectDemande.size()];
        vector < list<I_diarc *> > _vectLEntiers[vectDemande.size()];
        vector <int> _vectDemandeWLEntiers[vectDemande.size()];

//        _vectAreteHs.push_back(vector <simpleEdge> ());
//        for(int i=1; i<vectDemande.size(); i++)
//        {
//            _vectAreteHs.push_back(vector <simpleEdge> ());
//        }

//        #pragma omp parallel for private(dg)
        for(int i=0; i<vectDemande.size(); i++)
        {
            dg = grapheDemande[i];

            for(int r=0; r<(int)labelNumArcFracH[i].size()-1; r++)
            {
                for(int j=r+1; j<(int)labelNumArcFracH[i].size(); j++)
                {
                    list<I_dinode *> W2;
                    list<I_diarc *> lEntier;

                    if(dg->arcLier2(labelNumArcFracH[i][r],labelNumArcFracH[i][j],k,W2,lEntier))
                    {
                        //cout << "A1 = " << *arcFrac[r] << " et A2 = " << *arcFrac[j] << " ok" << endl;

                        simpleEdge se;
                        se.node1 = labelArcFracH[i*parDG+labelNumArcFracH[i][r]->num-1];
                        se.node2 = labelArcFracH[i*parDG+labelNumArcFracH[i][j]->num-1];
                        se.cap = 1;

                        _vectAreteHs[i].push_back(se);
                        _vectWHs[i].push_back(W2);
                        _vectLEntiers[i].push_back(lEntier);
                        _vectDemandeWLEntiers[i].push_back(i+1);
                    }
                    /*else
                    	//cout << "A1 = " << *arcFrac[r] << " et A2 = " << *arcFrac[j] << " non" << endl;*/
                }
            }
        }

        for(int i=0; i<vectDemande.size(); i++)
            for (int ii=0; ii<_vectAreteHs[i].size(); ii++)
            {
                vectAreteH.push_back(_vectAreteHs[i][ii]);
                vectWH.push_back(_vectWHs[i][ii]);
                vectLEntier.push_back(_vectLEntiers[i][ii]);
                vectDemandeWLEntier.push_back(_vectDemandeWLEntiers[i][ii]);
            }
    }

    /*I_Graph *gH = new I_Graph;

    InitGraph_from_list(vectAreteH,nbArcFracTotal,gH);
    DrawGraphViz(gH,"Test.dot",n,labelNumArcFracH2,nbDem);

    DeleteGraph(gH);
    delete gH;

    //sleep(10000);*/


	if(vectAreteH.size() > 0)
	{
        /**********************************/
		/*char s[256];
        static int qa = 0;

        qa++;
        cout << "QA == " << qa << endl;
        sprintf(s,".%d",qa);
        string str = "Test";
        str += s;
        str += ".dot";

		I_Graph *gH = new I_Graph;

		InitGraph_from_list(vectAreteH,nbArcFrac,gH);
		DrawGraphViz(gH,(char *)str.c_str(),n,
		labelNumArcFracH2,nbDem);

		DeleteGraph(gH);
		delete gH;*/
		/**********************************/

		//Vecteur de degre des sommets en terme d'arete noire
		//dans le graphe annexe. On ne s'interessera par la suite
		//qu'au sommet de degre non nul. On evitera ainsi de tomber
		//sur des cycles ayant deux aretes consecutives rouges
		vector<double> vectDegArcFrac(nbArcFracTotal);

		for(int i=0;i<nbArcFracTotal;i++)
			vectDegArcFrac[i] = 0;

		for(int i=0;i<vectAreteH.size();i++)
		{
			if(vectAreteH[i].cap == 1)
			{
				vectDegArcFrac[vectAreteH[i].node1-1]++;
				vectDegArcFrac[vectAreteH[i].node2-1]++;
			}
		}

		//Ce vecteur servira a marque les sommets qu'on trouve dans un cycle
		//On evitera de les prendre car en les reprenant on devrait retomber
		//sur le meme cycle en calculant un pcc.
		vector<bool> vectMarcArcFrac(nbArcFracTotal);

		for(int i=0;i<nbArcFracTotal;i++)
			vectMarcArcFrac[i] = false;

		//Calcul du graphe biparti
		//On construit une nouvelle liste d'aretes
		vector<simpleEdge> vectAreteBiparti;
		int n_bip = nbArcFracTotal;
		int BIG = 100000;

		for(int i=0;i<vectAreteH.size();i++)
		{
			simpleEdge se;

			if(vectAreteH[i].cap == 1)
			{
				se.node1 = vectAreteH[i].node1;
				se.node2 = vectAreteH[i].node2+n_bip;
				se.cap = 1;

				vectAreteBiparti.push_back(se);

				se.node1 = vectAreteH[i].node2;
				se.node2 = vectAreteH[i].node1+n_bip;
				se.cap = 1;

				vectAreteBiparti.push_back(se);
			}
			else
			{
				se.node1 = vectAreteH[i].node1;
				se.node2 = vectAreteH[i].node2;
				se.cap = BIG;

				vectAreteBiparti.push_back(se);

				se.node1 = vectAreteH[i].node1+n_bip;
				se.node2 = vectAreteH[i].node2+n_bip;
				se.cap = BIG;

				vectAreteBiparti.push_back(se);
			}
		}

		I_Graph *gBip = new I_Graph;

		InitGraph_from_list(vectAreteBiparti,2*nbArcFracTotal,gBip);

		//DrawGraphViz(gBip,"Test-Bip.dot",n);

		//On lance le calcul des pcc entre
		//chaque paire de sommet (u,u') du graphe biparti
//		#pragma omp parallel for
		for(int i=0;i<nbArcFracTotal;i++)
		{
			if(vectDegArcFrac[i] >= 1 /*&& !vectMarcArcFrac[i]*/)
			{

                double *T = NULL;
                int *P = NULL;
                int *PredArete = NULL;
                double *L = NULL;

                Dijkstra(gBip,i+1,L,P,PredArete,ALL_EDGES);

                //S'il existe un chemin de u->u' alors
                //il contient forcement un nombre impair
                //d'arete noire
                if(L[i+nbArcFracTotal] < INFINI /*&& (qa != 11 || (qa == 11 && i <= 7))*/)
                {
				    //cout << "Cycle Trouve avec I = " << i+1 << endl;

		            div_t chemin = div((int)L[i+nbArcFracTotal],BIG);

		            //cout << "NbImpair = " << chemin.rem << " Nb X = " << chemin.quot << endl;

		            //Construction de la contrainte correspondante
		            int Yu;
		            I_diarc *a;
		            list<I_diarc *> Ly1,Ly0;	//Liste des arcs de coeff 1 et 0
		            list<I_Edge *> Lx;	    //Liste des aretes dans la ctr
		            list<I_diarc*> F;			//Ensemble des arcs a enlever
		            vector< list<I_dinode *> > LWCtr;		//Liste des ensembles W utilise
		                                                //pour construire la contrainte

		            vector< list<I_diarc *> > LEntierCtr;	//Pour stocker les arcs entiers
		                                                //pris dans une coupe

		            vector<int> LDemandeWLEntierCtr;    //Pour stocker les indices des demandes associ�s � LW et LEntier

		            //Sert a tester si le cycle est
		            //elementaire ou non
		            vector<int> *testElementaire = new vector<int>;
		            double valCourArete;	//Sert aussi a tester si deux aretes consecutives
		                                    //ne sont pas noires. Si c'est le cas, on arete
		                                    //car on retrouvera un cycle plus correct ailleurs

		            bool terminer = false;
		            int numBip,numA,numW;
		            int som;
		            som =  i+1+nbArcFracTotal;
		            bool ok = true;

		            do
		            {
		                if(som > nbArcFracTotal)
		                    Yu = som-nbArcFracTotal;
		                else
		                    Yu = som;

		                /*if(som > nbArcFrac)
		                    cout << Yu << "'" << "->";
		                else
		                    cout << Yu << "->";*/

		                testElementaire->push_back(Yu);

		                /*cout << "BGIP 1 = " << gBip->Edges[PredArete[som-1]-1].cap << " " <<
		                gBip->Edges[PredArete[P[som-1]-1]-1].cap << endl;

		                if((P[som-1] != -1) && (gBip->Edges[PredArete[som-1]-1].cap == BIG) && (gBip->Edges[PredArete[P[som-1]-1]-1].cap == BIG) )
		                {
		                    cout << "Deux aretes noires consecutives" << endl;1738
		                    ok = false;
		                }*/

		                som = P[som-1];

		            }while((P[som-1] != -1) && (ok));


		            /*som =  i+1+nbArcFrac;
		            while(P[som-1] != -1)
		            {
		                if(som > nbArcFrac)
		                    Yu = som-nbArcFrac;
		                else
		                    Yu = som;

		                if(som > nbArcFrac)
		                    cout << Yu << "'" << "->";
		                else
		                    cout << Yu << "->";

		                som = P[som-1];
		            }

		            cout << i+1 << endl;*/

		            bool elementaire = true;
		            for(int j=0;j<(int)testElementaire->size();j++)
		            {
		                for(int r=j+1;r<(int)testElementaire->size();r++)
		                {
		                    if((*testElementaire)[j] == (*testElementaire)[r])
		                        elementaire = false;
		                }
		            }

		            delete testElementaire;

		            //Si le cycle est elementaire, on le
		            //garde. S'il est non elementaire, on
		            //ne le prend pas car, on sait qu'on
		            //retombera de toutes facon sur le
		            //cycle elementaire qu'il contient.


		            if(!elementaire);
		                //cout << "Cycle non elementaire" << endl;
		            else
		            {
		                //cout << "Cycle elementaire" << endl;

		                //On marque les sommets du cycle
		                som =  i+1+nbArcFracTotal;
		                while(P[som-1] != -1)
		                {
		                    if(som > nbArcFracTotal)
		                        Yu = som-nbArcFracTotal;
		                    else
		                        Yu = som;

		                    vectMarcArcFrac[Yu-1] = true;                                         ///Issue here probably?!

		                    som = P[som-1];
		                }

		                //On a un traitement si le premier arc qu'on rencontre
		                //est a 1 ou si le premier arc est a 0
		                som =  i+1+nbArcFracTotal;
		                numBip = PredArete[som-1];

		                //cout << "CAP = " << gBip->Edges[numBip-1].cap << endl;

		                if(gBip->Edges[numBip-1].cap == 1)
		                {
		                    do
		                    {
		                        numBip = PredArete[som-1];

		                        if(gBip->Edges[numBip-1].back->adjac->id != som)
		                            Yu = gBip->Edges[numBip-1].back->adjac->id;
		                        else
		                            Yu = gBip->Edges[numBip-1].adjac->id;

		                        if(Yu > nbArcFracTotal)
		                            Yu -= nbArcFracTotal;

		                        a = labelNumArcFracH3[Yu-1];
		                        Ly1.push_back(a);

		                        //S'il s'agit d'une arete rouge
		                        //alors, on supprime les deux derniers
		                        //sommets qui ont ete inseres
		                        //Et on ajoute leur arc correspondant
		                        //dans Lx. On les met eux meme dans Ly0
		                        if(gBip->Edges[numBip-1].cap != 1)
		                        {
		                            a = Ly1.back();
		                            Ly1.pop_back();
		                            Ly0.push_back(a);

		                            a = Ly1.back();
		                            Ly1.pop_back();
		                            Ly0.push_back(a);

		                            Lx.push_back(a->pedge);
		                        }

		                        //Si l'arete est noire, alors on
		                        //ajoute le W correspondant dans LW
		                        /*cout << "Cap 2 = " << gBip->Edges[numBip-1].cap << " " << numBip << endl;
		                        cout << gBip->Edges[numBip-1].back->adjac->id << " " <<
		                        gBip->Edges[numBip-1].adjac->id << endl; */

		                        if(gBip->Edges[numBip-1].cap == 1)
		                        {
		                            //numA = numBip/2;	//On prend la partie superieure
		                            numA = sup_part(numBip,2);

		                            numW = numA - pRouge;

		                            //Calcul du W correspondant
		                            LWCtr.push_back(vectWH[numW-1]);
		                            LEntierCtr.push_back(vectLEntier[numW-1]);
		                            LDemandeWLEntierCtr.push_back(vectDemandeWLEntier[numW-1]);
		                        }

		                        som = P[som-1];

		                    }while(P[som-1] != -1);
		                }
		                else
		                {
		                    //Dans ce cas, on commence par inserer
		                    //le premier sommet du cycle
		                    Yu = som;

		                    if(Yu > nbArcFracTotal)
		                        Yu -= nbArcFracTotal;

		                    a = labelNumArcFracH3[Yu-1];
		                    Ly1.push_back(a);

		                    do
		                    {
		                        numBip = PredArete[som-1];

		                        if(gBip->Edges[numBip-1].back->adjac->id != som)
		                            Yu = gBip->Edges[numBip-1].back->adjac->id;
		                        else
		                            Yu = gBip->Edges[numBip-1].adjac->id;

		                        if(Yu > nbArcFracTotal)
		                            Yu -= nbArcFracTotal;

		                        a = labelNumArcFracH3[Yu-1];
		                        Ly1.push_back(a);

		                        //S'il s'agit d'une arete rouge
		                        //alors, on supprime les deux derniers
		                        //sommets qui ont ete inseres
		                        //Et on ajoute leur arc correspondant
		                        //dans Lx. On les met eux meme dans Ly0
		                        if(gBip->Edges[numBip-1].cap != 1)
		                        {
		                            a = Ly1.back();
		                            Ly1.pop_back();
		                            Ly0.push_back(a);

		                            a = Ly1.back();
		                            Ly1.pop_back();
		                            Ly0.push_back(a);

		                            Lx.push_back(a->pedge);
		                        }

		                        //Si l'arete est noire, alors on
		                        //ajoute le W correspondant dans LW
		                        if(gBip->Edges[numBip-1].cap == 1)
		                        {
		                            numA = sup_part(numBip,2);	//On prend la partie superieure

		                            numW = numA - pRouge;

		                            LWCtr.push_back(vectWH[numW-1]);
		                            LEntierCtr.push_back(vectLEntier[numW-1]);
		                            LDemandeWLEntierCtr.push_back(vectDemandeWLEntier[numW-1]);
		                        }

		                        som = P[som-1];

		                        //Il faut s'arreter dans ce cas la
		                        //avant de revenir sur le dernier

		                    }while(P[som-1] != i+1);

		                    //Comme on s'est arrive avant la derniere arete du cycle
		                    //on doit ajouter le W de cette derniere. Bien sur, s'il s'agit
		                    //d'une arete noire
		                    //cout << "Som = " << som << endl;
		                    numBip = PredArete[som-1];

		                    if(gBip->Edges[numBip-1].cap == 1)
		                    {
		                        numA = sup_part(numBip,2);	//On prend la partie superieure
		                        numW = numA - pRouge;

		                        LWCtr.push_back(vectWH[numW-1]);
		                        LEntierCtr.push_back(vectLEntier[numW-1]);
		                        LDemandeWLEntierCtr.push_back(vectDemandeWLEntier[numW-1]);
		                    }
		                }

		                double poids_ca = 0;

		                //cout << "Ly1 = " << endl;
		                for(list<I_diarc *>::iterator it=Ly1.begin();it!=Ly1.end();it++)
		                {
		                    poids_ca += (*it)->cap;
		                    //cout << *(*it) << " " << endl;
		                }

		                //cout << "Ly0 = " << endl;
		                for(list<I_diarc *>::iterator it=Ly0.begin();it!=Ly0.end();it++)
		                {
		                    //cout << *(*it) << " " << endl;
		                }

		                //cout << "Lx = " << endl;
		                for(list<I_Edge *>::iterator it=Lx.begin();it!=Lx.end();it++)
		                {
		                    poids_ca += (*it)->X;
		                    //cout << "(" << (*it)->back->adjac->id << " " << (*it)->adjac->id << " " << (*it)->X << ")" << endl;
		                }

						int P_sz = LWCtr.size();
						int rhs_ca = sup_part(P_sz,2);

		                //cout << "POIDS CTR = " << poids_ca << " RHS = " << rhs_ca << endl;

		                //Arrive ici, on construit la contrainte
		                //Pour cela, il faut retrouver l'ensemble W
		                //correspondant a chaque arete noire

		                /********************************************************/
		                /* Ici on remplace les y restant par le X correspondant */

		                /*for(list<I_diarc *>::iterator it=Ly1.begin();it!=Ly1.end();it++)
		                {
		                    //cout << *(*it) << " " << endl;
		                    a = *it;
		                    Lx.push_back(a->pedge);
		                    Ly0.push_back(a);
		                }
		                Ly1.clear();*/
		                /**************************************************/

		                if(poids_ca < rhs_ca-EPSILON)
		                {
		                    /*for(list<I_diarc *>::iterator it=Ly1.begin();it!=Ly1.end();it++)
		                    {
		                        cout << *(*it) << " " << endl;
		                    }

		                    cout << "Ly0 = " << endl;
		                    for(list<I_diarc *>::iterator it=Ly0.begin();it!=Ly0.end();it++)
		                    {
		                        cout << *(*it) << " " << endl;
		                    }

		                    cout << "Lx = " << endl;
		                    for(list<I_Edge *>::iterator it=Lx.begin();it!=Lx.end();it++)
		                    {
		                        cout << "(" << (*it)->back->adjac->id << " " << (*it)->adjac->id << " " << (*it)->X << ")" << endl;
		                    }*/

                            /////////////////////////////////////////////////////////////////////

                            vector <int> *_cons;
                            double aCoeff = 0;

                            _cons = new vector <int> ();

                            _cons->push_back(0);
                            _cons->push_back(AG);
                            (*_cons)[0] += 2;

                            COUPE_AGGREGEE *new_ag;
                            new_ag = new COUPE_AGGREGEE(n,rhs_ca,&LWCtr,&LEntierCtr,&LDemandeWLEntierCtr,&Lx,&Ly0,&Ly1);

                            for(int i=0; i < gr->m_Edges; i++)
                            {

                                aCoeff = (int) new_ag->coeff(*c_info, i);

                                if (aCoeff > 0)
                                {
                                    _cons->push_back(aCoeff);
                                    _cons->push_back(i);
                                    (*_cons)[0] += 2;
                                }
                                //aCtrExpr += aCoeff * c_info->xSol[i];
                            }

                            _cons->push_back(GET);
                            _cons->push_back(rhs_ca);
                            (*_cons)[0] += 2;

                            seperated_cons.push_back( & ( (*_cons)[0] ) );
                            int n_p = 1;

                            res_ca++;

                            /////////////////////////////////////////////////////////////////////


		                    //Ajout de la contrainte dans le pool
//		                    ABA_BUFFER<ABA_CONSTRAINT *> p_constraints(master_,1);
//		                    COUPE_AGGREGEE *new_ca;
//
//		                    new_ca = new COUPE_AGGREGEE(getMaster(),n,rhs_ca,&LWCtr,&LEntierCtr,&LDemandeWLEntierCtr,&Lx,&Ly0,&Ly1);
//		                    p_constraints.push(new_ca);
//
//                            l1.lock();
//		                    addCons(p_constraints,getMaster()->getPoolCoupeAggregee(),0,0);
//
//		                    res_ca++;
//		                    l1.unlock();
		                }
		            }
		        }

		        delete[] T;
		        delete[] P;
		        delete[] PredArete;
		        delete[] L;

			}

            //cout << endl << endl;
        }

        DeleteGraph(gBip);
        delete gBip;


        //cout << "Nb CA = " << res_ca << endl;

        //sleep(36000);

        /*char c;
        cin >> c;*/
    }

    delete[] labelNumArcFracH;
    delete[] labelNumArcFracH2;


    return res_ca;
}



int extra_cuts::AjouteSPPartition(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie)
{
	bool sp_ok = false;
	Bool sp_part_result = false;  /// instead of defined False
	long nb_sp_part = 0;
	long n,m;
	I_Graph *gr_dem;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	gr_dem = &G_dem;

	long nb_chaine = 0,nb_chaine2 = 0;
	b_Node **chaine,**chaine2,*verif_part;
	long *chaine_sz,*chaine_sz2;
	long **sp_part_list;
	long *rhs;
	long *sp_p;

	chaine = new b_Node *[(n/2)];
	chaine_sz = new long[(n/2)];
	chaine2 = new b_Node *[(n/2)];
	chaine_sz2 = new long[(n/2)];
	rhs = new long[(n/2)];
	sp_p = new long[(n/2)];
	sp_part_list = new long *[(n/2)];


	for(int i=0;i<(n/2);i++)
	{
		chaine[i] = NIL_BN;
		chaine2[i] = NIL_BN;
		chaine_sz[i] = 0;
		chaine_sz2[i] = 0;

		rhs[i] = 0;
		sp_part_list[i] = NULL;
		sp_p[i] = 0;
	}

	//On calcule toutes les chaines qu'il y a dans le graphe
	get_sp_partition_chaine_steiner(gr,chaine,chaine_sz,&nb_chaine);

	//On construit des sous-chaines de longueur 2 et on verifie
	//si la SP-partition induite est violee ou non.
	nb_chaine2 = 0;
	for(int i=0;i<nb_chaine;i++)
	{
		b_Node *bncour = chaine[i];

		/*cout << "QQ == ";
		Print_b_Node_Set(chaine[i],stdout);*/

		while((bncour != NIL_BN) && (bncour->next != NIL_BN))
		{
			//On verifie que la partition induite est
			//admissible. Pour cela, comme les 2 premiers
			//elements de la partition sont reduits
			//chacun a un seul sommet, les sous-ensembles
			//induit sont d'office admissibles.
			//Il suffit donc de verifier que le troisieme
			//element: V\{V1,V2} est admissible

			verif_part = NIL_BN;
			I_Node *nptr = gr->Nodes[0].n_sh;

			do
			{
				if((nptr->id != bncour->id) && (nptr->id != bncour->next->id))
				{
					Insere_Node(&verif_part,nptr->id);
				}

				nptr = nptr->next_gr;

			}while(nptr != gr->Nodes[0].n_sh);

			//cin >> n;

			Bool admiss = admissible(gr_dem,verif_part);

			Delete_b_Node_Set(&verif_part);

			if(admiss)
			{
				//cout << "Partition admissible" << endl;

				Insere_Node(&chaine2[nb_chaine2],bncour->id);
				Insere_Node(&chaine2[nb_chaine2],bncour->next->id);

				chaine_sz2[nb_chaine2] = 2;

				/*cout << "Chaine " << i << endl;
				Print_b_Node_Set(chaine2[nb_chaine2],stdout);
				cout << endl;*/

				nb_chaine2++;
			}
			else
			{
				cout << "Non admissible" << endl;
			}

			bncour = bncour->next;
		}
	}

	/*********************************************************/
	//Destruction de la liste des chaines precedemment utilisee
	for(int i=0;i<nb_chaine;i++)
	{
		Delete_b_Node_Set(&chaine[i]);
	}


	delete[] chaine;
	delete[] chaine_sz;
	/**********************************************************/

	//sp_part_result = separation_sp_partition_2(gr,k,chaine,chaine_sz,nb_chaine,&sp_part_list,sp_p,&nb_sp_part,rhs);
	sp_part_result = separation_sp_partition_2(gr,k,chaine2,chaine_sz2,nb_chaine2,&sp_part_list,sp_p,&nb_sp_part,rhs);

	sp_ok = true;

	if(sp_part_result == true)   /// instead of defined True
	{
		/*Ajout de la contrainte dans le pool*/

		/////////////////////////////////////////////////////////////////////
        vector <int> *_cons;
        double aCoeff = 0;

        _cons = new vector <int> ();

        _cons->push_back(0);
        _cons->push_back(SPP);
        (*_cons)[0] += 2;

        SP_PARTITION *new_sp;
        int n_sp = 0;

        for(long i=0;i < nb_sp_part;i++)
		{
            new_sp = new SP_PARTITION(sp_part_list[i],sp_p[i],rhs[i]);

            for(int i=0; i < gr->m_Edges; i++)
            {

                aCoeff = (int) new_sp->coeff(*c_info, i);

                if (aCoeff > 0)
                {
                    _cons->push_back(aCoeff);
                    _cons->push_back(i);
                    (*_cons)[0] += 2;
                }
                //aCtrExpr += aCoeff * c_info->xSol[i];
            }

            _cons->push_back(GET);
            _cons->push_back(rhs[i]);
            (*_cons)[0] += 2;

            seperated_cons.push_back( & ( (*_cons)[0] ) );
            n_sp++;
		}

		/////////////////////////////////////////////////////////////////////

//		ABA_BUFFER<ABA_CONSTRAINT *> sp_constraints(getMaster(),nb_sp_part);
//		SP_PARTITION *new_sp;
//
//		for(long i=0;i<nb_sp_part;i++)
//		{
//			new_sp = new SP_PARTITION(getMaster(),sp_part_list[i],sp_p[i],rhs[i]);
//			sp_constraints.push(new_sp);
//		}
//
//		int n_sp = addCons(sp_constraints,getMaster()->getPoolSPPartition(),0,0);

		//addCons(sp_constraints);

		if(n_sp != nb_sp_part)
		{
			sp_ok = false;
		}
		else
		{
			sp_ok = true;
		}
	}
	else
	{
		nb_sp_part = 0;
	}


	for(int i=0;i<nb_chaine2;i++)
	{
		Delete_b_Node_Set(&chaine2[i]);
	}


	delete[] chaine2;
	delete[] chaine_sz2;
	delete[] sp_part_list;
	delete[] sp_p;
	delete[] rhs;

	if(sp_part_result == true)  /// instead of defined True
	{

		if(sp_ok)
		{
			//getMaster()->add_nb_sp_part(nb_sp_part);
			//cout << "vraieeeee"<< endl;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		nb_sp_part = 0;
	}

	return nb_sp_part;
}



int extra_cuts::AjouteSPPartition2(deque < int * > &seperated_cons, I_Graph *gr,int k,b_Node *frac_cycle,long cycle_sz,FILE *sortie)
{
	bool sp_ok = false;
	Bool sp_part_result = false ;  /// instead of defines False
	long nb_sp_part = 0;
	long n,m;
	I_Graph *gr_dem;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	gr_dem = &G_dem;

	long nb_chaine = 0,nb_chaine2 = 0;
	b_Node **chaine,**chaine2,*verif_part;
	long *chaine_sz,*chaine_sz2;
	long **sp_part_list;
	long *rhs;
	long *sp_p;

	chaine = new b_Node *[(n/2)];
	chaine_sz = new long[(n/2)];
	chaine2 = new b_Node *[(n/2)];
	chaine_sz2 = new long[(n/2)];
	rhs = new long[(n/2)];
	sp_p = new long[(n/2)];
	sp_part_list = new long *[(n/2)];


	for(int i=0;i<(n/2);i++)
	{
		chaine[i] = NIL_BN;
		chaine2[i] = NIL_BN;
		chaine_sz[i] = 0;
		chaine_sz2[i] = 0;

		rhs[i] = 0;
		sp_part_list[i] = NULL;
		sp_p[i] = 0;
	}

	//On calcule toutes les chaines qu'il y a dans le graphe
	get_sp_partition_chaine_steiner(gr,chaine,chaine_sz,&nb_chaine);

	//Si le graphe de demande est connexe, on cherche sur toute les chaines
	sp_part_result = false;  /// instead of defined False
	if(demandGraphConnected)
		sp_part_result = separation_sp_partition_2(gr,k,chaine,chaine_sz,nb_chaine,&sp_part_list,sp_p,&nb_sp_part,rhs);

	//Si le graphe de demandes n'est pas connexe ou si la précédente séparation
	//échoue alors on segment les en chaines de longueurs 2
	if(sp_part_result == false) /// instead of defined False
	{
		//On construit des sous-chaines de longueur 2 et on verifie
		//si la SP-partition induite est violee ou non.
		nb_chaine2 = 0;
		for(int i=0;i<nb_chaine;i++)
		{
			if(chaine_sz[i] >= 3)
			{
				b_Node *bncour = chaine[i];

				while((bncour != NIL_BN) && (bncour->next != NIL_BN))
				{
					Insere_Node(&chaine2[nb_chaine2],bncour->id);
					Insere_Node(&chaine2[nb_chaine2],bncour->next->id);

					chaine_sz2[nb_chaine2] = 2;
					nb_chaine2++;

					bncour = bncour->next;
				}
			}
		}

		sp_part_result = separation_sp_partition_2(gr,k,chaine2,chaine_sz2,nb_chaine2,&sp_part_list,sp_p,&nb_sp_part,rhs);
	}

	sp_ok = true;

	if(sp_part_result == true)  /// instead of defined True
	{
		/*Ajout de la contrainte dans le pool*/

		/////////////////////////////////////////////////////////////////////
        vector <int> *_cons;
        double aCoeff = 0;

        _cons = new vector <int> ();

        _cons->push_back(0);
        _cons->push_back(SPP);
        (*_cons)[0] += 2;

        SP_PARTITION *new_sp;
        int n_sp = 0;

        for(long i=0;i < nb_sp_part;i++)
		{
            new_sp = new SP_PARTITION(sp_part_list[i],sp_p[i],rhs[i]);

            for(int i=0; i < gr->m_Edges; i++)
            {

                aCoeff = (int) new_sp->coeff(*c_info, i);

                if (aCoeff > 0)
                {
                    _cons->push_back(aCoeff);
                    _cons->push_back(i);
                    (*_cons)[0] += 2;
                }
                //aCtrExpr += aCoeff * c_info->xSol[i];
            }

            _cons->push_back(GET);
            _cons->push_back(rhs[i]);
            (*_cons)[0] += 2;

            seperated_cons.push_back( & ( (*_cons)[0] ) );
            n_sp++;
		}

		/////////////////////////////////////////////////////////////////////

//		ABA_BUFFER<ABA_CONSTRAINT *> sp_constraints(getMaster(),nb_sp_part);
//		SP_PARTITION *new_sp;
//
//		for(long i=0;i<nb_sp_part;i++)
//		{
//			new_sp = new SP_PARTITION(getMaster(),sp_part_list[i],sp_p[i],rhs[i]);
//			sp_constraints.push(new_sp);
//		}
//
//		int n_sp = addCons(sp_constraints,getMaster()->getPoolSPPartition(),0,0);
		//addCons(sp_constraints);

		if(n_sp != nb_sp_part)
		{
			sp_ok = false;
		}
		else
		{
			sp_ok = true;
		}
	}
	else
	{
		nb_sp_part = 0;
	}

	for(int i=0;i<nb_chaine;i++)
	{
		Delete_b_Node_Set(&chaine[i]);
	}

	for(int i=0;i<nb_chaine2;i++)
	{
		Delete_b_Node_Set(&chaine2[i]);
	}

	delete[] chaine;
	delete[] chaine_sz;
	delete[] chaine2;
	delete[] chaine_sz2;
	delete[] sp_part_list;
	delete[] sp_p;
	delete[] rhs;

	if(sp_part_result == true)  /// instead of defined True
	{
		if(sp_ok)
		{
			//getMaster()->add_nb_sp_part(nb_sp_part);
			cout << "SP2 vraieeeee"<< endl;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		nb_sp_part = 0;
	}

	return nb_sp_part;
}

