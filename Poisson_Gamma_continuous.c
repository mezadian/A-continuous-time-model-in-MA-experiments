#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <time.h>
#define MAXTYPES 1501   // maximum number of fitness classes to keep track of
#define MAXGENS  30     // maximum number of generations of growth
#define MAXMUTES 10001  // maximum number of entries in the DFE

int main(int argc, char** argv) {

    srand(time(NULL));
    clock_t begin = clock();
    FILE* fpout, * fps, * fpfit;
    char filename[100];
    long Nvisible;   // minimum number of indls in a visible colony
    int i, j, ngens, iline, nlines, itransfer, ntransfers, nbabies, T;
    int initialanc, ianc, nNtypes;
    int Nend, nlucky, ilucky, Ntmp;
    double Wsave, dfe[3][MAXMUTES], fmax, gammax, gamshift, delw;
    double r, r_1, r_2 = 0.0, r_3 = 0.0;
    int gamk, focal_row, bin_index=0, mu_row, mu_index=0, jsave = 0;
    double mutotal;
    double selective_effect[MAXTYPES - 1], death_mu[MAXTYPES - 1], rates[MAXTYPES - 1];
    double A[MAXTYPES - 1], cumlSum[MAXTYPES - 1], newA[MAXTYPES - 1];
    long Ntypes_total[MAXTYPES - 1], N;
    double t = 0.0, death = 0.2, sumA, time2nextevent;
    int nsinglets = 0;   // total number of single-mutation entries in DFE
    
    long seed = -32;
    if (argc > 1) seed = -((long)(atof(argv[1])));
    float ran1(long*), gamdev(int, long*), poidev(float, long*);

    long** Ntypes = (long**)malloc(2 * sizeof(long*));
    for (i = 0;i <= 1;i++) Ntypes[i] = (long*)malloc(MAXTYPES * sizeof(long));

    Nvisible = 1e5;
    mutotal = 0.043;
    nNtypes = MAXTYPES-1;
    ngens = 15;
    ntransfers = 5;
    nlines = 1000;
    //T = 0; // average final colony size for ngens

    fprintf(stdout, "%d generations between transfers\n", ngens);

    fmax = (double)(1.0 / (1.0 - death)) - 1.0;
    //printf("fmax:%f\n", fmax);
    delw = (double)((1.0 + fmax) / nNtypes);
    //printf("delw:%f\n", delw);
    initialanc = (int)(1.0 / delw);    // initial ancestor has relative fitness 1
    //printf("initialanc: %d\n", initialanc);


    gamk = 10;   // shape parameter for gamma distn of underlying DFE
    gammax = 0.5;
    gamshift = (gammax - 0.25); // We can then shift the mean to below 1 by this amount

    for (int i = 0; i < nNtypes; i++) {
        selective_effect[i] = -1.0 + i * delw;
        death_mu[i] = 1.0 - (1.0 + selective_effect[i]) * (1.0 - death);
        rates[i] = 1.0;
        //printf("mutotal[1200]: %f\n", mutotal);
        //printf("death_mu: %f\n");
        //printf("selective_effect, death_mu, rates: %f %f %f\n", selective_effect[i],death_mu[i], rates[i]);
    }
    
    //printf("death_mu[ianc], selective_effect[ianc]: %f %f\n", death_mu[1200], selective_effect[1200]);
    //printf("death_mu[nNtypes-1], selective_effect[nNtypes-1]: %f %f\n", death_mu[nNtypes-1], selective_effect[nNtypes-1]);

    // check that we can open the output file properly before doing the work
    sprintf_s(filename, 100, "madata/bd_%dtrans_%dlines_%dk_%1.2fM_%1.5fmu_%dgens_%dseed.out", ntransfers, nlines, gamk, gammax - gamshift, mutotal, ngens, -seed);
    int pout = fopen_s(&fpout, filename, "w");
    if (pout != 0) { fprintf(stdout, "Unable to open output file\n"); exit(1); }
    sprintf_s(filename, 100, "madata/bd_%dtrans_%dlines_%dk_%1.2fM_%1.5fmu_%dgens_%dseed_fitness.out", ntransfers, nlines, gamk, gammax - gamshift, mutotal, ngens, -seed);
    int pfit = fopen_s(&fpfit, filename, "w");
    if (pfit != 0) { fprintf(stdout, "Unable to open output file\n"); exit(1); }

    for (iline = 1;iline <= nlines;iline++) {
        //printf("iline: %d\n", iline);
        ianc = initialanc;

        for (itransfer = 1;itransfer <= ntransfers;itransfer++) {
            //printf("itransfer: %d\n", itransfer);

            //do {

                // initially start with one indl with ancestar fitness
                for (i = 0;i < nNtypes;i++) {
                    Ntypes[0][i] = 0; Ntypes[1][i] = 0;
                }
                
                Ntypes[0][ianc] = 1;
                //printf("Ntypes[0][ianc]: %d %d\n", Ntypes[0][ianc],ianc);
                Wsave = ianc * delw;   // save the ancestor fitness here
                //printf("Wsave: %f\n", Wsave);


                t = 0.0;

                while (t < ngens) {
                    for (int j = 0; j < nNtypes;j++) {
                        Ntypes_total[j] = Ntypes[0][j] + Ntypes[1][j];
                        //printf("Ntypes_total: %d\n", Ntypes_total[j]);
                    }
                    //printf("Ntypes_total[ianc]: %d\n", Ntypes_total[ianc]);
                    for (int j = 0; j < nNtypes; j++) {
                        A[j] = Ntypes_total[j] * rates[j];
                        //printf("A: %f\n", A[j]);

                        if (j == 0) {
                            cumlSum[j] = A[j];
                        }
                        else {
                            cumlSum[j] = cumlSum[j - 1] + A[j];
                        }
                        //printf("cumlSum: %f\n", cumlSum[j]);
                    }
                    //printf("A[ianc]: %f\n", A[ianc]);
                    //printf("cumlSum[ianc]: %f\n", cumlSum[ianc]);

                    time2nextevent = -log(ran1(&seed)) / cumlSum[nNtypes - 1];
                    //printf("time2nextevent: %f\n", time2nextevent);
                    //exit(1); 

                    for (int j = 0; j < nNtypes; j++) {
                        newA[j] = (double)cumlSum[j] / cumlSum[nNtypes - 1];
                        //printf("newA: %f\n", newA[j]);
                    }
                    //printf("newA[ianc]: %f\n", newA[ianc]);

                    r_1 = ran1(&seed);
                    bin_index = 0;
                    while (r_1 > newA[bin_index]) { bin_index++; }
                    //printf("r_1, bin_index: %f %d\n", r_1, bin_index);

                    // Generate r_3 and ensure it satisfies one of the conditions

                    r_3 = ran1(&seed) * (Ntypes[0][bin_index] + Ntypes[1][bin_index]);
                    if (r_3 < Ntypes[0][bin_index]) {
                        focal_row = 0; // Condition satisfied, assign 0
                    }
                    else {
                        focal_row = 1; // Condition satisfied, assign 1
                    }
                    //printf("r_3, focal_row: %f %d\n", r_3, focal_row);

                    // The expected # of wild-type offspring follow by fission distribution
                    // which is binomial distribution with N=2.
                    // nbabies gives us 0, 1, or 2 offspring
                    nbabies = poidev((float)(2 * (1 - death_mu[bin_index])), &seed);
                    //printf("nbabies: %d\n", nbabies);
                    
                    // the following line means that after having nbabies 
                    // the mother cell will eliminate
                    Ntypes[focal_row][bin_index] -= 1;

                    // Once the mother cell give birth and eliminate from the lineage 
                    // their offspring can mutant or not at that instant.

                    for (i = 0; i < nbabies; i++) {
                        r_2 = ran1(&seed);
                        //printf("r_2: %f\n", r_2);
                        if (r_2 < mutotal) {
                            // The next 2 lines tell us to put any new indls into row 1
                            // of nextNtypes, unless they are new mutations of the ancestor,
                            // in which case they will be stored in row 0
                            mu_index = 0;
                            mu_row = 1;
                            if (bin_index == ianc) mu_row = 0;

                            do {
                                // Recalculate r and mu_index
                                r = 1.0 + gammax - gamshift - gammax * gamdev(gamk, &seed) / (float)gamk;
                                //r =  gamdev(gamk, &seed);
                                //fprintf(fpfit, " %f\n", r);
                                //newW = bin_index*delw*r;  mu_index = newW/delw;
                                mu_index = (int)(bin_index * r);
                                //printf("mu_index, bin_index: %d %d\n", mu_index, bin_index);
                            } while (mu_index < 0 || mu_index >= nNtypes);

                            //the above do-while loop means that if mu_index is out of range 
                            // [0,nNtypes] do again until it satify

                            Ntypes[mu_row][mu_index] += 1;

                            //printf("mutation, mu_row, mu_index:%d %d %d\n", Ntypes[mu_row][mu_index], mu_row, mu_index);
                        }
                        else {
                            Ntypes[focal_row][bin_index] += 1;
                            //printf("birth, focal_row, bin_index:%d %d %d\n", Ntypes[focal_row][bin_index], focal_row, bin_index);
                        }
                    }

                    N = 0;
                    for (j = 0; j < nNtypes; j++) { N += Ntypes[0][j] + Ntypes[1][j]; }
                    //printf("N: %d\n", N);

                    if (N == 0) {
                        //printf("Extinction: all elements of Ntypes are zero, restart the loop.\n");
                        t = 0.0; // Reset the time
                        Ntypes[0][ianc] = 1; //reinitialize ancestor
                    }

                    else { t = t + time2nextevent; }
                    //printf("t: %f\n", t);

                } // loop on t, growth is finished

                // Nend gives total size of colony at end of growth
                Nend = 0;
                for (i = 0;i < nNtypes;i++)
                {
                    Nend += Ntypes[0][i] + Ntypes[1][i];
                }
                //printf("Nend: %d\n", Nend);
           // } while (Nend < T);

            //printf("Nend: %d\n", Nend);

            //if (ngens > 20 && Nend < Nvisible)
                //fprintf(stderr, "Warning: we chose a colony with only %d indls\n", Nend);
            nlucky = (long)(ran1(&seed) * Nend);
            //printf("nlucky: %d\n", nlucky);
            Ntmp = 0; jsave = 1;
            ilucky = ianc;
            for (i = 0;i < nNtypes;i++) {
                for (j = 0;j <= 1;j++) {
                    Ntmp += Ntypes[j][i];
                    if (Ntmp > nlucky) {
                        ilucky = i;
                        jsave = j;
                        i = nNtypes + 1; j = 2;  // jump out of both loops
                    }
                }
            }
            //printf("Ntmp: %d\n", Ntmp);
            //printf("ilucky, jsave: %d %d\n", ilucky, jsave);

            /* add this mutation to the list if previous ancestor
                       and this lucky differ by  a single mutation  */
                       //printf("ianc: %d\n", ianc);
            if ((ilucky != ianc) && (jsave == 0)) {
                nsinglets++;
                //printf("nsinglets: %d\n", nsinglets);
                if (nsinglets > MAXMUTES) fprintf(stderr, "Error: MAXMUTES exceeded\n");
                dfe[0][nsinglets] = ilucky * delw / Wsave;
                //printf("dfe[0][nsinglets]: %f %d\n", dfe[0][nsinglets], nsinglets);
                dfe[1][nsinglets] = (double)itransfer;
                //printf("dfe[1][nsinglets]: %f\n", dfe[1][nsinglets]);
                dfe[2][nsinglets] = ilucky * delw;
                //printf("dfe[2][nsinglets]: %f %d\n", dfe[2][nsinglets], nsinglets);
            }

            // now lucky becomes the new ancestor for the next growth
            ianc = (int)ilucky;
            //printf("new ancestor: % d\n", ianc);
            for (i = 0;i < nNtypes;i++) {
                Ntypes[0][i] = Ntypes[1][i] = 0;
            }
            Ntypes[0][ianc] = 1; // this type is the new ancestor
            Wsave = ianc * delw;

        } // loop on itransfers, do the next growth/sampling now

    } // at this point, this line is done, return to do a new line

    fprintf(fpfit, "\n");
    fprintf(stdout, "%d ", iline);

    // at this point, all the lines are done.  Save results.
    fprintf(stdout, "\n");
    for (i = 1;i <= nsinglets;i++)
    fprintf(fpout, "%f %f %f\n", dfe[0][i], dfe[1][i], dfe[2][i]);
    fclose(fpout);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    fprintf(stdout, "Time spent: %f\n", time_spent);
}