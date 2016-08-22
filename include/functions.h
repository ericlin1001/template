#ifdef FUNCTION_H
#else
#define FUNCTION_H
class Function{
#define MAX_FUNCTION_NAME 150
	private:
		char shortName[50];
		char funName[MAX_FUNCTION_NAME];
		double xlow,xup;
		double fbest;
		bool isFindMin;
		int numDim;
		//
		int feCounter;
	private:
	public:
		static double u(double x,double a,double k,double m){
			if(x>a)return k*pow(x-a,m);
			if(x<-a)return k*pow(-x-a,m);
			return 0;
		}
		virtual double operator()(const double *xs,int size){
			feCounter++;
			return 0;
		}
		inline double f(const vector<double>&xs){
			return operator()(&xs[0],xs.size());
		}
	public:
		Function(const char *s,double xlow,double xup,double fbest,bool isFindMin,int numDim){
			this->xlow=xlow;
			this->xup=xup;
			this->fbest=fbest;
			this->isFindMin=isFindMin;
			this->numDim=numDim;
			strcpy(shortName,s);
			if(isFindMin){
				sprintf(funName,"{%s(%f,%f)^%d fmin:%f}",s,xlow,xup,numDim,fbest);
			}else{
				sprintf(funName,"{%s(%f,%f)^%d fmax:%f}",s,xlow,xup,numDim,fbest);
			}
			feCounter=0;
		}
		int getfeCounter()const{return feCounter;}
		double getBest()const{return fbest;}
		bool getIsFindMin()const{return isFindMin;}
		inline bool isFBetter(double a,double b){
			return isFindMin^(a>=b);
		}
		int getNumDim()const{return numDim;}
		double getRange(int botOrUp){
			if(botOrUp==0)return xlow;
			return xup;
		}
		const char *getShortName()const{return shortName;}
		const char *getName()const{return funName;}
};

#define DefFunction(name,xlow,xup,fbest,isFindMin) class name : public Function{\
	public: name(int numDim):Function(#name,xlow,xup,fbest,isFindMin,numDim){}\
			virtual double operator()(const double *xs,int size){\
				Function::operator()(xs,size);
#define EndDef }};	
DefFunction(PDEF3,-10,10,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double t=xs[i+1]-xs[i];
		double t1=xs[i]-1;
		res+=t*t*100.0+t1*t1;
	}
	return res;
EndDef

DefFunction(PDEF4,-500,500,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		res+=-xs[i]*sin(sqrt(fabs(xs[i])));
	}
	res+=(double)418.9829*(double)size;
	return res;
EndDef

DefFunction(F1,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double x=xs[i];
		res+=x*x;
	}
return res;
	EndDef

DefFunction(F2,-10,10,0,true)
	double res=0.0;
	double sum=0.0;
	double mul=1.0;
	for(int i=0;i<size;i++){
		double fabsx=fabs(xs[i]);
		sum+=fabsx;
		mul*=fabsx;
	}
res=sum+mul;
return res;
EndDef

DefFunction(F3,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double insum=0.0;
		for(int j=0;j<=i;j++){
			insum+=xs[j];
		}
		res+=insum*insum;
	}
return res;
EndDef

DefFunction(F4,-100,100,0,true)
	double res=fabs(xs[0]);
	for(int i=1;i<size;i++){
		double tmp=fabs(xs[i]);
		if(tmp<res)res=tmp;
	}
return res;
EndDef

//untest:
DefFunction(F5,-30,30,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double tmp=pow(xs[i+1]-xs[i]*xs[i],2)*100.0+pow(xs[i]-1.0,2);
		res+=tmp;
	}
return res;
EndDef

DefFunction(F6,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		int tmp=floor(xs[i]+0.5);
		res+=tmp*tmp;
	}
return res;
EndDef

DefFunction(F7,-1.28,1.28,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],4)*(double)(i+1);
		res+=tmp;
	}
	res+=drand();
return res;
EndDef

DefFunction(F8,-500,500,-12569.5,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=-xs[i]*sin(sqrt(fabs(xs[i])));
		res+=tmp;
	}
return res;
EndDef

DefFunction(F9,-5.12,5.12,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],2)-(double)10.0*cos(xs[i]*2.0*M_PI)+10.0;
		res+=tmp;
	}
return res;
EndDef

DefFunction(F10,-32,32,0,true)
	double res=0.0;
	double sumx2=0.0;
	double sumcosx=0.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		sumcosx+=cos(xs[i]*M_PI*2.0);
	}
res=-20.0*exp(-0.2*sqrt(sumx2/(double)size))-exp(sumcosx/(double)size)+20.0+M_E;
return res;
EndDef

DefFunction(F11,-600.0,600.0,0,true)
	double res=0.0;
	double sumx2=0.0;
	double mulcos=1.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		mulcos*=cos(xs[i]/sqrt((double)i+1));
	}
res=sumx2/4000.0-mulcos+1.0;
return res;
EndDef

DefFunction(F12,-50,50,0,true)
	double res=0.0;
	double y1=1.0+(xs[0]+1.0)/4.0;
	double yd=1.0+(xs[size-1]+1.0)/4.0;
	double sumy=0.0;
	double sumu=0.0;
	//
	double yi,yi1;
	yi=y1;
	for(int i=0;i<size-1;i++){
		yi1=1.0+(xs[i+1]+1.0)/4.0;
		sumy+=pow(yi-1.0,2)*(1.0+10.0*pow(sin(M_PI*yi1),2));
		yi=yi1;
	}
for(int i=0;i<size;i++){
	sumu+=Function::u(xs[i],10,100,4);
}
res=M_PI/(double)size*(10.0*pow(sin(M_PI*y1),2)+sumy+pow(yd-1,2))
	+sumu;
	return res;
	EndDef

DefFunction(F13,-50,50,0,true)
	double res=0.0;
	double sumx=0.0;
	double sumu=0.0;
	for(int i=0;i<size-1;i++){
		sumx+=pow(xs[i]-1,2)*(1+pow(sin(3.0*M_PI*xs[i+1]),2));
	}
for(int i=0;i<size;i++){
	sumu+=Function::u(xs[i],5,100,4);
}
double xd=xs[size-1];
res=0.1*(pow(sin(3.0*M_PI*xs[0]),2)+sumx+
		pow(xd-1.0,2)*(1+pow(sin(2.0*M_PI*xd),2)))+sumu;
return res;
EndDef
class FunctionFactory{
	private:
		vector<Function*>fs;
		FunctionFactory(int numDim){
			fs.resize(4);
			fs[0]=new F1(numDim);
			fs[1]=new F3(numDim);
			fs[2]=new PDEF3(numDim);
			fs[3]=new PDEF4(numDim);
		}
		static FunctionFactory*instance;
	public:
		static FunctionFactory &Instance(int numDim){
			if(instance==0)instance=new FunctionFactory(numDim);
			return *instance;
		}
		/*
		   void setNumDim(int numDim){
		   }
		 */
		Function*getFunction(int index)const{
			return fs[index];
		}
		int getNumFunction()const{
			return fs.size();
		}
		~FunctionFactory(){
			for(int i=0;i<getNumFunction();i++){
				delete fs[i];
			}
		}
};
FunctionFactory*FunctionFactory::instance=0;
#endif
