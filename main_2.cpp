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
	if(isnan(expo)){
		cout<<u<<" "<<v<<" "<<theta<<endl;
		exit(0);
	}
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
	if(num<=0){
		cout<<"PDF<0   ";
		cout<<u<<" "<<v<<" "<<theta<<" "<<(num/(4*eps*eps))<<endl;
		exit(0);
	}
	return (num/(4*eps*eps));
}


//function for recursively multiplying Bivariate copulas given
// constant v and array u,theta
// TODO: Add check for copula families
double recursiveMultiply(double u[assets],double theta1[assets], double theta2[assets],double v1, double v2, char family='g'){
	double mult=1.0;
	double eval=0.0;
	for(int i=0;i<assets;i++){
		if(family=='g'){
		//	cout<<"Calculating Gumbel"<<endl;
			double inner=gumbel(u[i],v1,theta1[i]);
		//	cout<<u[i]<<" "<<theta1[i]<<" "<<v1<<" "<<inner<<endl;
			eval=gumbelPdf(inner,v2,theta2[i]);
		}
		else if(family=='j'){
			double inner=joe(u[i],v1,theta1[i]);
			eval=joePdf(inner,v2,theta2[i]);
		}
		mult*=eval;
	}
	
	if(isnan(mult)){
		cout<<v1<<" "<<v2<<" "<<endl;
		exit(0);
	}
	
	return mult;
}

double gaussLegendre2D(double u[assets], double theta1[assets], double theta2[assets], char family='g'){
	double sum=0.0;
	for(int i=0;i<integPoints;i++){
		for(int j=0;j<integPoints;j++){
			sum=sum+gaussLegendWeights[i]*gaussLegendWeights[j]*recursiveMultiply(u,theta1,theta2,gaussLegendPoints[j],gaussLegendPoints[i],family);
		}
	}
	return sum;	
}

double logLikelihood(double data[points][assets], double theta1[assets], double theta2[assets],char family='g'){
	double sum=0.0;
	for(int i=0;i<points;i++){
		double g=gaussLegendre2D(data[i],theta1,theta2,family);
		if(isnan(g)||g<=0){
			//cout<<g<<endl;
			cout<<"Point "<<i+1<<" : "<<g<<endl;
			exit(0);
		}
		sum+=log(g);
	}
	return sum;
}



int iterations=50;
double alpha=0.001; //Learning Rate

void find_optimal_parameters(double[points][assets],double theta1[assets],double theta2[assets]){
	
	double thetaUpd1[assets];
	double thetaFinal1[assets];
	double thetaUpd2[assets];
	double thetaFinal2[assets];
	for(int i=0;i<iterations;i++){
		double l=logLikelihood(U,theta1,theta2);
		cout<<"Iteration "<<i+1<<"    Likelihood = "<<l<<"   theta1 = {";
		for(int j=0;j<assets;j++){
			thetaUpd1[j]=theta1[j];
			cout<<theta1[j]<<",";
		}
		cout<<"}    theta2 = {";
		for(int j=0;j<assets;j++){
			thetaUpd2[j]=theta2[j];
			cout<<theta2[j]<<",";
		}
		cout<<"}"<<endl;
		
		for(int j=0;j<assets;j++){
			
			thetaUpd1[j]+=eps;
			double derivative=(logLikelihood(U,thetaUpd1, theta2)-l)/eps;
			thetaUpd1[j]-=eps;
			double change=alpha*derivative;
			if(change<0){
				thetaFinal1[j]=theta1[j]+max(-5.0,change);
			}
			else{
				thetaFinal1[j]=theta1[j]+min(5.0,change);
			}
			if(thetaFinal1[j]<=1.0){
				thetaFinal1[j]=1.0001;
			}
		}
		for(int j=0;j<assets;j++){
			
			thetaUpd2[j]+=eps;
			double derivative=(logLikelihood(U,theta1,thetaUpd2)-l)/eps;
			thetaUpd2[j]-=eps;
			double change=alpha*derivative;
			if(change<0){
				thetaFinal2[j]=theta2[j]+max(-5.0,change);
			}
			else{
				thetaFinal2[j]=theta2[j]+min(5.0,change);
			}
			if(thetaFinal2[j]<=1.0){
				thetaFinal2[j]=1.0001;
			}
		}
		theta1=thetaFinal1;
		theta2=thetaFinal2;
	}	
}

int main(){
	read_data();
	read_GLData();
	//cout<<U[0][0]<<endl;
	double theta[assets];
	for(int i=0;i<assets;i++){
		theta[i]=3;
	}
	find_optimal_parameters(U,theta,theta);
	
	
//	for(int i=0;i<integPoints;i++){
//		cout<<gaussLegendPoints[i]<<" "<<gaussLegendWeights[i]<<endl;
//	}
	//cout<<logLikelihood(U,theta,theta);
	
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
