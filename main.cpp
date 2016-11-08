#include<bits/stdc++.h>
#include<fstream>

using namespace std;

const int assets=15;
const int points=250;
double U[points][assets];
double eps=0.00001;
const int integPoints=17; // Gauss Legendre points
double gaussLegendPoints[integPoints];
double gaussLegendWeights[integPoints];

void read_data(){
	ifstream ifs("snp500.csv");
	string st;
	for(int i=0;i<points;i++){
		getline(ifs,st);
	    istringstream ss1(st);
	    for(int j=0;j<assets;j++){
	    	ss1>>U[i][j];
	    }
	}
}

void read_GLData(){
	ifstream ifs("gldata");
	string st;
	for(int i=0;i<integPoints;i++){
		getline(ifs,st);
	    istringstream ss1(st);
	    ss1>>gaussLegendPoints[i]>>gaussLegendWeights[i];
	}
}

double gumbel(double u,double v, double theta){
	double first=pow(-log(u),theta);
	double second=pow(-log(v),theta);
	double expo=pow(first+second,(1/theta));
	return exp(-expo);
}

double joe(double u,double v, double theta){
	double first=pow(1-u,theta);
	double second=pow(1-v,theta);
	double expo=pow(first+second-first*second,(1/theta));
	return 1-expo;
}

double joePdf(double u, double v, double theta){
	double umin,umax,vmin,vmax;
	umin=u-eps;
	vmin=v-eps;
	umax=u+eps;
	vmax=v+eps;
	if(u-eps<0.0){
		umin=0.0;
	}
	if(u+eps>1.0){
		umax=1.0;
	}
	if(v-eps<0.0){
		vmin=0.0;
	}
	if(v+eps>1.0){
		vmax=1.0;
	}
	double num=joe(umax,vmax,theta)-joe(umax,vmin,theta)-joe(umin,vmax,theta)+joe(umin,vmin,theta);
	return (num/(4*eps*eps));
}

double gumbelPdf(double u, double v, double theta){
	double umin,umax,vmin,vmax;
	umin=u-eps;
	vmin=v-eps;
	umax=u+eps;
	vmax=v+eps;
	if(u-eps<0.0){
		umin=0.0;
	}
	if(u+eps>1.0){
		umax=1.0;
	}
	if(v-eps<0.0){
		vmin=0.0;
	}
	if(v+eps>1.0){
		vmax=1.0;
	}
	double num=gumbel(umax,vmax,theta)-gumbel(umax,vmin,theta)-gumbel(umin,vmax,theta)+gumbel(umin,vmin,theta);
	if(num<=0){
		return 1e-7;
	}
	return (num/(4*eps*eps));
}


//function for recursively multiplying Bivariate copulas given
// constant v and array u,theta
// TODO: Add check for copula families
double recursiveMultiply(double u[assets],double theta[assets],double v, char family='g'){
	double mult=1.0;
	double eval=0.0;
	for(int i=0;i<assets;i++){
		if(family=='g'){
			eval=gumbelPdf(u[i],v,theta[i]);
		}
		else if(family=='j'){
			eval=joePdf(u[i],v,theta[i]);
		}
		mult*=eval;
	}
	return mult;
}

double gaussLegendre1D(double u[assets], double theta[assets], char family='g'){
	double sum=0.0;
	for(int i=0;i<integPoints;i++){
		sum=sum+gaussLegendWeights[i]*recursiveMultiply(u,theta,gaussLegendPoints[i],family);
	}
	return sum;	
}

double logLikelihood(double data[points][assets], double theta[assets], char family='j'){
	double sum=0.0;
	for(int i=0;i<points;i++){
		double g=gaussLegendre1D(data[i],theta,family);
		//cout<<"Point "<<i+1<<" : "<<g<<endl;
		sum+=log(g);
	}
	return sum;
}



int iterations=300;
double alpha=0.0001; //Learning Rate

void find_optimal_parameters(double[points][assets],double theta[assets]){
	
	double thetaUpd[assets];
	double thetaFinal[assets];
	for(int i=0;i<iterations;i++){
		double l=logLikelihood(U,theta);
		cout<<"Iteration "<<i+1<<"    Likelihood = "<<l<<"   theta = {";
		for(int j=0;j<assets;j++){
			thetaUpd[j]=theta[j];
			cout<<theta[j]<<",";
		}
		cout<<"}"<<endl;
		
		for(int j=0;j<assets;j++){
			
			thetaUpd[j]+=eps;
			double derivative=(logLikelihood(U,thetaUpd)-l)/eps;
			thetaUpd[j]-=eps;
			double change=alpha*derivative;
			if(change<-1.0){
				change=-1.0;
			}
			else if(change>1.0){
				change=1.0;
			}
			thetaFinal[j]=theta[j]+change;
			if(thetaFinal[j]<=1.0){
				thetaFinal[j]=1.0001;
			}
		}
		theta=thetaFinal;
	}	
}

int main(){
	read_data();
	read_GLData();
	double theta[assets];
	for(int i=0;i<assets;i++){
		theta[i]=3;
	}
	find_optimal_parameters(U,theta);
	
//	for(int i=0;i<integPoints;i++){
//		cout<<gaussLegendPoints[i]<<" "<<gaussLegendWeights[i]<<endl;
//	}
	//cout<<logLikelihood(U,theta);
	
	/*  
	Test for gumbel PDF
	-------------------------
	double theta=1.71;
	cout<<"start"<<endl;
	for(int i=1;i<100;i++){
		double u,v;
		u=0.01*i;
		for(int j=1;j<100;j++){
			v=0.01*j;
			double res=gumbelPdf(u,v,theta);
			if(res<0){
				cout<<"U ="<<u<<"     V ="<<v<<"       pdf ="<<res<<endl;
			}
		}
	
	--------------------------
	*/
	
	
	return 0;
}
