#

conformityDecomposition = function (adj, Cl = NULL) 
{
      if (  is.null(dim(adj) )) stop("Input adj is not a matrix or data frame. ")
    if ( dim(adj)[[1]] < 3)   stop("The adjacency matrix has fewer than 3 rows. This network is trivial and will not be evaluated.")
  if (!.is.adjmat(adj)) 
        stop("The input matrix is not a valid adjacency matrix!")
    diag(adj) = 0
    if (!is.null(Cl)) {
        if (length(Cl) != dim(adj)[[1]]) {
            stop(paste("The length of the class assignment Cl does not equal the number", 
                "of rows of the adjcency matrix. length(Cl) != dim(adj)[[1]]. \n", 
                "Something is wrong with your input"))
        }
     if (sum(is.na(Cl))>0 ) stop("Cl must not contain missing values (NA)." ) 
    }

A.CF=matrix(0, nrow=dim(adj)[[1]], ncol= dim(adj)[[2]] )
   diag(A.CF)=1
   if ( is.null(Cl) )  {
    Conformity = .NPC.iterate(adj)$v1
    if (sum(adj^2,na.rm=T)==0) {Factorizability=NA} else {
 A.CF=outer(Conformity, Conformity) - diag(Conformity^2)
    Factorizability = 1 - sum((adj - A.CF)^2)/sum(adj^2)} 
diag(A.CF)=1
    output = list(A.CF=A.CF, Conformity=data.frame( Conformity),  IntermodularAdjacency=1, Factorizability = Factorizability)
   } 

   if ( !is.null(Cl) )  {
   Cl=factor(Cl)
   Cl.level=levels( Cl )
   if ( length(Cl.level)>100 ) warning(paste("Your class assignment variable Cl contains",  length(Cl.level), "different classes. I assume this is a proper class assignment variable. But if not, stop the calculation, e.g. by using the Esc key on your keybord."))
   
   Conformity=rep(NA, length(Cl) )
   listConformity=list()
  IntramodularFactorizability=rep(NA, length(Cl.level) )
IntermodularAdjacency= matrix(0,nrow=length(Cl.level),ncol=length(Cl.level) )
diag(IntermodularAdjacency)=1
IntermodularAdjacency=data.frame(IntermodularAdjacency)
dimnames(IntermodularAdjacency)[[1]]=as.character(Cl.level)
dimnames(IntermodularAdjacency)[[2]]=as.character(Cl.level)

numeratorFactorizability=0
   for (i in 1:length(Cl.level) ) {
    restclass= Cl== Cl.level[i]
    if (sum(restclass)==1) { A.help=0; CF.help =0;   Conformity[restclass]=CF.help; 
A.CF[restclass,restclass]=CF.help*CF.help - CF.help^2
 }
if (sum(restclass)==2) {
A.help=adj[restclass,restclass];diag(A.help)=0 
CFvalue=sqrt(adj[restclass,restclass][1,2]);
CF.help= c(CFvalue , CFvalue )
Conformity[restclass]=CF.help
A.CF[restclass,restclass]=outer(CF.help, CF.help) - diag(CF.help^2)
 }

  if (sum(restclass)>2) {
  A.help=adj[restclass,restclass];diag(A.help)=0 ;
CF.help = .NPC.iterate(A.help )$v1 
Conformity[restclass]=CF.help
A.CF[restclass,restclass]=outer(CF.help, CF.help) - diag(CF.help^2)
  } 
if (length(CF.help)>1) {numeratorFactorizability= numeratorFactorizability+sum(   (A.help-outer(CF.help,CF.help)+diag(CF.help^2) )^2 )}
   
 listConformity[[i]] =CF.help 
      if (sum(A.help^2,na.rm=T)==0 |  length(CF.help)==1 ) {IntramodularFactorizability[i]=NA} else {
    IntramodularFactorizability[i] = 1 - sum((A.help - outer(CF.help, CF.help) + 
        diag(CF.help^2))^2)/sum(A.help^2,na.rm=T)  }
   } # end of for loop over i
  
   if ( length(Cl.level)==1) {IntermodularAdjacency[1,1]=1} else {
   for (i in 1:(length(Cl.level)-1)  ) {
   for (j in (i+1):length(Cl.level) ) {
   restclass1= Cl== Cl.level[i]
   restclass2= Cl== Cl.level[j]
   A.inter=adj[restclass1,restclass2]  
   mean.CF1=mean(listConformity[[i]], na.rm=T)
   mean.CF2=mean(listConformity[[j]], na.rm=T)
   if (   mean.CF1* mean.CF2    != 0 ) {
   IntermodularAdjacency[i,j]= mean(A.inter,na.rm=T)/(mean.CF1* mean.CF2)    
    IntermodularAdjacency[j,i]= IntermodularAdjacency[i,j]
      }

 if (  length(listConformity[[i]])==1 |  length(listConformity[[j]])==1    )   {
numeratorFactorizability=  numeratorFactorizability+ 2*sum( (A.inter- IntermodularAdjacency[i,j]* listConformity[[i]] * listConformity[[j]])^2  ) 
A.CF[restclass1,restclass2]= IntermodularAdjacency[i,j]* listConformity[[i]] *listConformity[[j]] 
A.CF[restclass2,restclass1]= IntermodularAdjacency[j,i]* listConformity[[j]] *listConformity[[i]]   
  } else {
numeratorFactorizability=  numeratorFactorizability+ 2*sum( (A.inter- IntermodularAdjacency[i,j]*outer( listConformity[[i]] , listConformity[[j]]) )^2  )
A.CF[restclass1,restclass2]= IntermodularAdjacency[i,j]* outer(listConformity[[i]], listConformity[[j]]  )
A.CF[restclass2,restclass1]= IntermodularAdjacency[j,i]* outer(listConformity[[j]], listConformity[[i]]  )

 } # end of else
   } # end of for (i in
   } # end of for (j in
   diag(adj)=NA
 Factorizability=  1-numeratorFactorizability/ sum(adj^2,na.rm=T)  

   } # end of if else statement
   
diag(A.CF)=1
    output = list(A.CF=A.CF,  Conformity=Conformity, IntermodularAdjacency= IntermodularAdjacency, Factorizability = Factorizability,  Cl.level =Cl.level, IntramodularFactorizability=  IntramodularFactorizability, listConformity=listConformity )
   } # end of   if ( !is.null(Cl) ) 
    output
  }
