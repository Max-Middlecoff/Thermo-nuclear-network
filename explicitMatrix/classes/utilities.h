
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>

using namespace std;
using std::string;

/* Class Utilities to hold utility useful utility functions.  Functions are
 * declared static so that they can be invoked without having to instantiate
 * objects of type Utilities.  For example, Utilities::returnNetIndexZN (Z, N).
 * We will often have other classes inherit from Utilities so that they can
 * access its static functions.  This class definition must precede definitions of 
 * classes that inherit from it.
 */

class Utilities{
    
    private:
    
    public:
        
        // Static function Utilities::showTime() to return date and local time as a 
        // character array
        
        static char* showTime(){
            time_t now = time(0);         // Current date/time
            char* tnow = ctime(&now);     // convert to string
            return tnow;
        }
        
        // -------------------------------------------------------------------------
        // Static function Utilities::interpolate_T(t) to find an interpolated T9
        // as a function of time if hydroProfile = true.
        // -------------------------------------------------------------------------
        
        static double interpolate_T(double t, double T9_start){
            
            // Will call spline interpolator in hydro profile table to return 
            // T9 at this value of time t.  For now, we just return T9_start.
            
            double temp;            // Interpolated temperature
            temp = T9_start;        // Temporary constant value for testing
            
            // **********
            // Call to spline interpolator for temperature goes here
            // **********
            
            return temp;
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::interpolate_rho(t) to find an interpolated 
        // density rho as a function of time if hydroProfile is false.
        // -------------------------------------------------------------------------
        
        static double interpolate_rho(double t, double rho_start){
            
            // Will call spline interpolator in hydro profile table to return 
            // rho at this value of time t.  For now, we just return rho_start.
            
            double rhonow;             // Interpolated density
            rhonow = rho_start;        // Temporary constant value for testing
            
            // **********
            // Call to spline interpolator for density goes here
            // **********
            
            return rhonow;
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::log10Spacing() to find num equal log10 
        // spacings between two numbers start and stop. The equally log-spaced
        // numbers are placed in the array v passed with the pointer *v.
        // -------------------------------------------------------------------------
        
        static void log10Spacing(double start, double stop, int num, double *v){
            
            double logtmin = log10(start);
            double logtmax = log10(stop);
            double tempsum = logtmin;
            double expofac = (logtmax - logtmin) / (double) num;
            v[0] = pow(10, logtmin);
            for(int i=0; i<num; i++){
                tempsum += expofac;
                v[i] = pow(10, tempsum);
            }
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::plotOutput() to output data at regular 
        // intervals of the integration to a file suitable for plotting. Assumes
        // the existence of a subdirectory gnu_out. May crash if this directory
        // does not exist. Assuming gnuplot for plotting, but the output file
        // is whitespace-delimited ascii, so any plotting program could be used
        // to read it. Lines beginning with # are comments in gnuplot.
        // -------------------------------------------------------------------------
        
        static void plotOutput( 
                                FILE * pFileD, 
                                bool doASY, 
                                bool plotFluxes, 
                                bool doPE, 
                                bool dopf, 
                                int totalTimeSteps, 
                                int totalTimeStepsZero, 
                                int plotSteps,
                                double tplot[],
                                double dtplot[],
                                double slowestRatePlot[],
                                char reacLabel[][35], // @todo: SHOULD NOT DO THIS!!!!!
                                int slowestRateIndexPlot[], 
                                double fastestRatePlot[],
                                int fastestRateIndexPlot[],
                                double dt_FEplot[],
                                double dt_EAplot[],
                                double dt_trial[],
                                char isoLabel[][5],  // @todo: SHOULD NOT DO THIS!!!!!
                                double EReleasePlot[],
                                double dEReleasePlot[],
                                int numAsyplot[],
                                int ISOTOPES,
                                int numRG_PEplot[],
                                int numberRG,
                                double sumXplot[],
                                double Xplot[][200], // @todo: SHOULD NOT DO THIS!!!!!
                                double FplusSumPlot[][200], // @todo: SHOULD NOT DO THIS!!!!!
                                double FminusSumPlot[][200] // @todo: SHOULD NOT DO THIS!!!!!
                              ) { 

            // Open files for ascii output. Assumes that the subdirectory
            // gnu_out already exists. If it doesn't, will compile but
            // may crash when executed.
            
            FILE * pFile;
            pFile = fopen("gnu_out/gnufile.data","w");
            
            FILE * pFile2;
            pFile2 = fopen("gnu_out/gnufile2.data","w");
            
            FILE * pFile3;
            pFile3 = fopen("gnu_out/gnufileFlux.data","w");
            
            // Following array controls which mass fractions are exported to plotting
            // file.  The entries in plotXlist[] are the species indices for the
            // isotopes in the network to be plotted. For small networks export all;
            // for large networks we will usually export only a representative subset.
            // Hardwired for now, but eventually we should read the entries of this
            // array in from a data file.
            
            //             int maxPlotIsotopes = 16;
            //             int plotXlist[maxPlotIsotopes];
            //             for(int i=0; i<maxPlotIsotopes; i++){
            //                 plotXlist[i] = i;
            //             }
            

            //             int plotXlist[] = {0,1,2,3,4,5,6};                              // pp
               int plotXlist[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};      // alpha
            //             int plotXlist[] = {0,1,2,3};                                    // 4-alpha
            //             int plotXlist[] = {0,1,2};                                      // 3-alpha
            //             int plotXlist[] = {0,1,2,3,4,5,6,7};                            // cno
            //             int plotXlist[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};      // cnoAll
            //             
            //             int plotXlist[] = 
            //             {4,12,20,28,35,42,52,62,72,88,101,114,128,143,0,1,
            //             13,16,43,49,147,132,123,38,25,32,30,34,7,18,21,38};   // 150-isotope select)
            
            
            // Get length LX of array plotXlist holding the species indices for
            // isotopes that we will plot mass fraction X for.
            
            int LX = sizeof(plotXlist)/sizeof(plotXlist[0]);
            
            string str1 = "#    t     dt     |E|  |dE/dt| Asy  Equil  sumX";
            string strflux = "\n#    t     dt   ";
            string app = "  ";
            string app1;
            string appflux;
            string Xstring = "X(";
            string Fpstring = "F+(";
            string Fmstring = "F-(";
            string dFstring = "dF(";
            string iso;
            
            fprintf(pFileD, "\n\n");
            
            if(doASY){
                fprintf(pFile, "# ASY");
                fprintf(pFile2, "# ASY");
                if(plotFluxes){fprintf(pFile3, "# ASY");}
                fprintf(pFileD, "# ASY");
            } else {
                fprintf(pFile, "# QSS");
                fprintf(pFile2, "# QSS");
                if(plotFluxes){fprintf(pFile3, "# QSS");}
                fprintf(pFileD, "# QSS");
            }
            
            if(doPE){
                fprintf(pFile, "+PE");
                fprintf(pFile2, "+PE");
                if(plotFluxes){fprintf(pFile3, "+PE");}
                fprintf(pFileD, "+PE");
            } 
            
            if(dopf){
                fprintf(pFile, " method (with partition functions): ");
                fprintf(pFile2, " method (with partition functions): ");
                if(plotFluxes){fprintf(pFile3, " method (with partition functions): ");}
                fprintf(pFileD, "+ method (with partition functions): ");
            } else {
                fprintf(pFile, " method (no partition functions): ");
                fprintf(pFile2, " method (no partition functions): ");
                if(plotFluxes){fprintf(pFile3, " method (no partition functions): ");}
                fprintf(pFileD, "+ method (no partition functions): "); 
            }
            
            fprintf(pFile, "%d integration steps ", totalTimeSteps - totalTimeStepsZero);
            fprintf(pFile2, "%d integration steps ", totalTimeSteps - totalTimeStepsZero);
            if(plotFluxes) fprintf(pFile3, "%d integration steps ", 
                totalTimeSteps - totalTimeStepsZero);
            fprintf(pFileD, "%d integration steps ", totalTimeSteps - totalTimeStepsZero);
                
            FPRINTF_CPU;
            FPRINTF_CPU2;
            FPRINTF_CPUD;
            
            fprintf(pFile, "# All quantities except Asy, RG_PE, and sumX are log10(x)\n");
            fprintf(pFile, "# Log of absolute values for E and dE/dt as they can be negative\n");
            fprintf(pFile, "# Units: t and dt in s; E in erg; dE/dt in erg/g/s; others dimensionless \n");
            fprintf(pFile, "#\n");
            
            string str2 = "#      t       dt  2/Rmin   Reaction_Rmin    1/Rmax   Reaction_Rmax";
            str2 += ("     dt_FE   dt_EA   trial_dt\n");
            fprintf(pFile2, "# All double quantities are log10(x); rates in units of s^-1\n#\n");
            fprintf(pFile2, stringToChar(str2));
            
            for(int i=0; i<plotSteps; i++){
                
                fprintf(pFile2, "%7.4f %7.4f %7.4f %s %7.4f %s %7.4f %7.4f %7.4f\n", 
                    tplot[i], dtplot[i], log10(1.0/slowestRatePlot[i]), 
                    reacLabel[ slowestRateIndexPlot[i]],
                    log10(2.0/fastestRatePlot[i]),
                    reacLabel[ fastestRateIndexPlot[i]],
                    log10(dt_FEplot[i]), log10(dt_EAplot[i]), log10(dt_trial[i])
                );
            }
            
            // Write header for file pointed to by pFile
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[plotXlist[i]];
                app.append(Xstring);
                app.append(iso);
                app.append(")     ");
            }
            
            str1.append(app);
            str1.append("\n");
            fprintf(pFile, stringToChar(str1));
            
            fprintf(pFile, "\n");
            
            // Write header for file pointed to by pFile3

            for(int i=0; i<LX; i++){
                iso = isoLabel[plotXlist[i]];
                appflux.append(Fpstring);
                appflux.append(iso);
                appflux.append(")   ");
            }
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[plotXlist[i]];
                appflux.append(Fmstring);
                appflux.append(iso);
                appflux.append(")   ");
            }
            
            for(int i=0; i<LX; i++){
                iso = isoLabel[plotXlist[i]];
                appflux.append(dFstring);
                appflux.append(iso);
                appflux.append(")   ");
            }
            
            strflux.append(appflux);
            fprintf(pFile3, stringToChar(strflux));
            
            
            // Loop over timesteps for plot output writing the data to the file 
            // line by line using concatenated fprintf statements.
            
            for(int i=0; i<plotSteps; i++){
                
                // Initial data fields for t, dt, sumX, fraction of asymptotic
                // isotopes, and fraction of reaction groups in equilibrium.
                
                fprintf(pFile, "%+6.3f %+6.3f %6.3f %6.3f %5.3f %5.3f %5.3f",
                    tplot[i], dtplot[i], EReleasePlot[i], dEReleasePlot[i], 
                    (double)numAsyplot[i]/(double)ISOTOPES,
                    (double)numRG_PEplot[i]/(double)numberRG,
                    sumXplot[i]
                );
                
                // Now add one data field for each X(i) in plotXlist[]. Add
                // 1e-24 to X in case it is identically zero since we are
                // taking the log.
                
                for(int j=0; j<LX; j++){
                    fprintf(pFile, " %5.3e", log10(Xplot[plotXlist[j]][i]+1e-24));
                }
                
                fprintf(pFile, "\n");
                
                // Fluxes
                
                fprintf(pFile3, "\n%+6.3f %+6.3f", tplot[i], dtplot[i]);
                
                // Now add one data field for each FplusSumPlot. Add
                // 1e-24 to X in case it is identically zero since we are
                // taking the log.
                
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", log10(abs( FplusSumPlot[j][i]+1e-24) ));
                }
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", log10(abs( FminusSumPlot[j][i]+1e-24)));
                }
                for(int j=0; j<LX; j++){
                    fprintf(pFile3, " %5.3e", 
                        log10( abs(FplusSumPlot[j][i] - FminusSumPlot[j][i] + 1e-24) ));
                }
                
                cout.flush();   // Force buffer dump
                
            }
            
            // Close output files
            
            fclose (pFile);
            fclose (pFile2);
            fclose (pFile3);
            
        }   // End plotOutput()
        
        
        // Static function Utilities::plotHydroprofile() to send hydro profile 
        // to plotting file. Only invoked if hydroProfile and plotHydroProfile 
        // are true.
        
        static void plotHydroProfile(
                                        int hydroLines,
                                        double hydroTime[],
                                        double hydroTemp[],
                                        double hydroRho[]
                                    ){
            
            FILE * pHydro;
            pHydro = fopen("gnu_out/hydroProfile.out","w");
            fprintf(pHydro, "#  time          T          rho");
            
            for (int i=0; i<hydroLines; i++){
                fprintf(pHydro, "\n%6.4e  %6.4e  %6.4e", 
                    hydroTime[i], hydroTemp[i], hydroRho[i]);
            }
            
            fclose (pHydro);
        }
        
        
        // -------------------------------------------------------------------------
        // Static function Utilities::sumMassFractions() to return the current sum of the
        // mass fractions X(i) in the network. If the network conserves particle
        // number this sum should be equal to 1.0.
        // -------------------------------------------------------------------------
        
        static double sumMassFractions( int isotopes, double X[] ) {
            
            double sum = 0.0;
            for(int i=0; i<isotopes; i++){
                sum += X[i];
            }
            return sum; 
        }
        
        
        // ------------------------------------------------------------------
        // Static function Utilities::sumXEquil() to return the sum of mass 
        // fractions for isotopes participating in at least one RG currently
        // in partial equilibrium
        // ------------------------------------------------------------------
        
        static double sumXEquil( int isotopes, double X[], bool isotopeInEquil[] ) {
            
            double sum = 0;
            for(int i=0; i<isotopes; i++){
                if( isotopeInEquil[i] ){
                    sum += X[i];
                }
            }
            return sum;
        }
        
        
        // -----------------------------------------------------------------------
        // Static function Utilities::sumXNotEquil() to return the sum of mass 
        // fractions for isotopes not participating in any RG currently in
        // partial equilibrium
        // -----------------------------------------------------------------------
        
        static double sumXNotEquil( int isotopes, double X[], bool isotopeInEquil[] ) {
            
            double sum = 0;
            for(int i=0; i<isotopes; i++){
                if( !isotopeInEquil[i] ){
                    sum += X[i];
                }
            }
            return sum;
        }
        
    
        // -------------------------------------------------------------------------
        // Static function Utilities::returnNetIndexZN(Z,N) to return the network 
        // vector index given Z and N for the isotope.  Returns -1 if no match.
        // -------------------------------------------------------------------------
        
        static int returnNetIndexZN(int z, int n, int numberSpecies, int Z[], int N[] ) {
            
            for (int i = 0; i < numberSpecies; i++) {
                if (Z[i] == z && N[i] == n) return i;
            }
            return -1;
        }
        
        
        // -----------------------------------------------------------------------
        // Static function Utilities::returnNetIndexSymbol(char* symbol) to 
        // return the network vector index given the symbol for the isotope.
        // -----------------------------------------------------------------------
        
        static int returnNetIndexSymbol(char* symbol, int numberSpecies, char isoLabel[][5]) { // @todo DO NOT DO THIS!!!!!
            
            int result;
            for (int i = 0; i < numberSpecies; i++) {
                result = strcmp(isoLabel[i], symbol);
                if (result == 0){
                    return i;
                }
            }
            return -1;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::isInNet(Z,N) to return true if given (Z,N) 
        // is in the network, false otherwise.
        // ----------------------------------------------------------------------
        
        static bool isInNet(int z, int n, int numberSpecies, int Z[], int N[]) {
            
            if (returnNetIndexZN(z, n, numberSpecies, Z, N) < 0) {
                return false;
            } else {
                return true;
            }
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::minimumOf(x,y) to return minimum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int minimumOf(int x, int y) {            // integers
            return (x < y) ? x : y;
        }
        
        static double minimumOf(double x, double y) {   // doubles
            return (x < y) ? x : y;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::maximumOf(x,y) to return maximum of two 
        // numbers.  Overloaded to accept either integer or double arguments.
        // ----------------------------------------------------------------------
        
        static int maximumOf(int x, int y) {            // integers
            return (x > y) ? x : y; 
        }
        
        static double maximumOf(double x, double y){    // doubles
            return (x > y) ? x : y; 
        }
        
        
        // ----------------------------------------------------------------------
        // Using the C++ class string instead of char to handle strings 
        // won't compile on my system unless #include <iostream> is included. The
        // C printf command also won't work correctly because it isn't typesafe.
        // See the discussion at
        //
        //    https://stackoverflow.com/questions/10865957/printf-with-stdstring
        //
        // Instead, print with the cout command (which requires #include <iostream>)
        // Example: 
        //
        //      string test("howdy");
        //      string test2(" do");
        //      cout << "Your string is " << test + test2;
        //
        // If instead you try
        //
        //      printf("\n%s", test+test2);
        //
        // it will likely print garbage. However, you can print a string with printf
        // if it is first converted to a char:
        //
        //      string s = "Howdy World!";
        //      char cs[s.size() + 1];
        //      strcpy(cs, &s[0]);	// or strcpy(cs, s.c_str());
        //      printf("\n\nstring=%s\n", strcpy(cs, &s[0]));
        //  
        // Typically a string type can be printed with cout but a string given 
        // to printf usually displays garbage because of type issues in the 
        // C function printf noted above. The function stringToChar(string) defined
        // below converts a string to a corresponding character array, which either 
        // printf or cout can print. Presently assumes the string has no more than 50 
        // characters.  Change the dimension of cs[] to change that.
        // ----------------------------------------------------------------------
        
        static char* stringToChar(string s){
            static char cs[50];
            strcpy(cs, &s[0]);     // alternatively strcpy(cs, s.c_str());
            return cs;
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::startTimer() to start a timer.  Stop timer
        // and display elapsed time with Utilities::stopTimer().
        // ----------------------------------------------------------------------
        
        static void startTimer(){
            START_CPU     // Start a timer for rate calculation
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::stopTimer() to stop timer started with 
        // Utilities::startTimer() and display results.
        // ----------------------------------------------------------------------
        
        static void stopTimer(){
            
            STOP_CPU;        // Stop the timer
            printf("\n");
            PRINT_CPU;       // Print timing information for rate calculation
            printf("\n");
        }
        
        
        // ----------------------------------------------------------------------
        // Static function Utilities::testTimerCPU() to test CPU timer by 
        // executing long, pointless loop.
        // ----------------------------------------------------------------------
        
        static void testTimerCPU(){
            
            double a, b;
            
            START_CPU;
            for (long count = 1l; count < 500000l; count++) {
                a = sqrt(count);
                b = 1.0/logf(a);
                a = logf(b)/sqrt(a);
            }
            STOP_CPU;
            PRINT_CPU_TEST;
            printf("\n");
        }
    
};    // End class Utilities
