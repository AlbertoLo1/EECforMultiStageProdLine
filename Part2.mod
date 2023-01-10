/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Alberto Loffredo
 * Creation Date: 8 nov 2021 at 11:19:43
 *********************************************/


int nStates = ...;

range States = 1..nStates;
range NrStations = 1..2;
float NrMachines[NrStations] =...;
float BuffCap[NrStations] =...;

float Lambda = ...;
float Alpha[States] = ...;

string pos1 = ...;
string pos2 = ...;
string pos3 = ...;
string pos4 = ...;
string pos5 = ...;
string pos6 = ...;

float OF;

float Target = 0;

float Costs[States][States] = ...;
float Prob[States][States] = ...;
float Mach[States][NrStations] = ...;
 
dvar float+ x[States][States];
 
maximize sum(j in States) sum(l in States)  x[j][l] * Costs[j][l];
 
subject to{

forall(m in States)
  sum(i in States) x[m][i] - sum(i in States) Lambda*Prob[i,m]*x[i][m] == Alpha[m];  

forall(jj in NrStations)
sum (m in States, i in States) x[m][i]*Mach[i][jj] >= Target * NrMachines[jj] * sum(m in States, i in States)x[m][i];

} 

execute {

OF = cplex.getObjValue();

}