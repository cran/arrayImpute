arrayImpute<-function()
{
   require(RGtk2)
   require(cairoDevice)
   require(impute)
   require(matlab)
   require(quantreg)
   require(gnm)
   require(nlme)
   K<<-10
   Result.Save<<-NULL
   high.list<<-NULL
   high.chip.list<<-NULL
   gene.location<<-NULL
   true.flag<<-FALSE

#############################################################################################
## control GeneList, ChipList, and output Table

geneListSelectRow.handler <- function(w, r, c, e, u = NULL)
{     high.list<<-sort(cbind(high.list,(r+1)))
      high.list<<-unique(high.list)
}

geneListUnselectRow.handler <- function(w, r, c, e, u = NULL)
{     high.list<<-high.list[high.list!=(r+1)]
      high.list<<-unique(high.list)
}

chipListSelectRow.handler <- function(w, r, c, e, u = NULL)
{    high.chip.list<<-sort(cbind(high.chip.list,(r+1)))
      high.chip.list<<-unique(high.chip.list)
}

chipListUnselectRow.handler <- function(w, r, c, e, u = NULL)
{     high.chip.list<<-high.chip.list[high.chip.list!=(r+1)]
      high.chip.list<<-unique(high.chip.list)
}

showMissingList<- function(missing.Name)
{    if(ncol(missing.Name)==3)
     {  temp.p<-9
        temp<-matrix(rep(" ",temp.p*nrow(missing.Name)),ncol=temp.p)
        missing.Name<-cbind(missing.Name,temp)
     }
     for(i in 1: nrow(missing.Name))
        gtkCListAppend(geneList,as.character(missing.Name[i,]))
}

showChipList<- function()
{   chip.name<<-as.matrix(colnames(data.missing))
     for(i in 1: nrow(chip.name))
        gtkCListAppend(ChipList,as.character(chip.name[i,]))
}

removeMissingList <- function(n)
{   for(i in (n-1):0)
       gtkCListRemove(geneList,i)
}

warning.gtk <- function(messages) 
{  warn.win <- gtkWindow(show = FALSE)
    box <- gtkVBox(TRUE, 0)
    box$SetUsize(-1, 80)
    box$PackStart(gtkLabel(messages))
    box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
    box.small$PackStart(button <- gtkButton("Ok"), expand = TRUE, fill = FALSE)
    button$SetUsize(60, 25)
    button$AddCallback("clicked", function(x) {
      warn.win$Hide()
      warn.win$Destroy()
    })
    warn.win$Add(box)
    warn.win$Show()
}
  
############################################################################################# 
### In File menu
## Read microaray data with missing

Mnu.read.file.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open Microarray data file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         data.missing<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
         rownames(data.missing)<<-as.character(data.missing[,1])
         data.missing<<-data.missing[,-1]
         missing.indicator(data.missing)
         AVE.impute<<-apply(data.missing[missing.I[,2],],1,function(x){mean(x,na.rm=T)})
         missing.Name<<-cbind(missing.Name,round(AVE.impute,3))
         showMissingList(missing.Name)
         showChipList()
         contents <- paste(" Read ", exp.info.name, sep=" ")
         contents <- paste(contents,"\n",sep=" ")
         Result.Save<<-c(Result.Save,contents)
         gtkTextBufferInsertAtCursor(new.text,contents)
         gtkTextViewSetBuffer(outputText,new.text)

         win$Hide()
         win$Destroy()
    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })
}

## Read Gene Location with pin information

Mnu.read.location.file.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open Gene Location file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         gene.location<<-NULL
         gene.location.o<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
	 pin<-paste(gene.location.o[,1],gene.location.o[,2],sep="")
         n<-table(pin)[1]
         p<-length(table(gene.location.o[,3]))
         q<-length(table(gene.location.o[,4]))
         x<-gene.location.o[,4]+q*(gene.location.o[,2]-1)
         y<-gene.location.o[,3]+p*(gene.location.o[,1]-1)
         gene.location<<-cbind(x,y)
	 x<<-x
	 y<<-y
         win$Hide()
         win$Destroy()

    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })
}

missing.indicator<-function(data.missing)
{  missing.I<<-NULL
   missing.Name<<-NULL
   gene.name<-rownames(data.missing)
   chip.name<-colnames(data.missing)
   for(k in 1:ncol(data.missing))
   {  temp.I<-which(is.na(data.missing[,k]))
      missing.Name<<-rbind(missing.Name,cbind(rep(chip.name[k],length(temp.I)),gene.name[temp.I]))
      missing.I<<-rbind(missing.I,cbind(rep(k,length(temp.I)),temp.I))
   }
   colnames(missing.Name)<<-c("Chip","Gene")
   colnames(missing.I)<<-c("Chip","Gene")
   missing.Name.temp<<-missing.Name
   miss.temp<<-missing.I
   kNN.impute<<-rep(" ",nrow(missing.Name.temp))
   LLS.impute<<-rep(" ",nrow(missing.Name.temp))
   RLSP.impute<<-rep(" ",nrow(missing.Name.temp))
   BPCA.impute<<-rep(" ",nrow(missing.Name.temp))
   User1.impute<<-rep(" ",nrow(missing.Name.temp))
   User2.impute<<-rep(" ",nrow(missing.Name.temp))
}

## Save Result

Mnu.save.result.handler<-function(w,u=NULL)
{  if(kNN.impute[1]!=" ")
   {   full.data<-data.missing
       full.data[is.na(data.missing)]<-kNN.impute
       win<- gtkFileSelection("Save kNN imputation Result as")
       ok <- win[["OkButton"]]
       cancel <- win[["CancelButton"]]
       ok$AddCallback("clicked", function(x) {
            file.name <- win$GetFilename()
            write.table(full.data,file.name,quote=FALSE,sep=',')
            win$Hide()
            win$Destroy()
       }) 
    }
   if(LLS.impute[1]!=" ")
   {   full.data<-data.missing
       full.data[is.na(data.missing)]<-LLS.impute
       win<- gtkFileSelection("Save LLS imputation Result as")
       ok <- win[["OkButton"]]
       cancel <- win[["CancelButton"]]
       ok$AddCallback("clicked", function(x) {
            file.name <- win$GetFilename()
            write.table(full.data,file.name,quote=FALSE,sep=',')
            win$Hide()
            win$Destroy()
       }) 
    }
   if(RLSP.impute[1]!=" ")
   {   full.data<-data.missing
       full.data[is.na(data.missing)]<-RLSP.impute
       win<- gtkFileSelection("Save RLSP imputation Result as")
       ok <- win[["OkButton"]]
       cancel <- win[["CancelButton"]]
       ok$AddCallback("clicked", function(x) {
            file.name <- win$GetFilename()
            write.table(full.data,file.name,quote=FALSE,sep=',')
            win$Hide()
            win$Destroy()
       }) 
    }

   if(BPCA.impute[1]!=" ")
   {   full.data<-data.missing
       full.data[is.na(data.missing)]<-BPCA.impute
       win<- gtkFileSelection("Save BPCA imputation Result as")
       ok <- win[["OkButton"]]
       cancel <- win[["CancelButton"]]
       ok$AddCallback("clicked", function(x) {
            file.name <- win$GetFilename()
            write.table(full.data,file.name,quote=FALSE,sep=',')
            win$Hide()
            win$Destroy()
       }) 
    }
}

## Exit

Mnu.exit.handler<-function(w,u=NULL)
{  main$Hide()
   main$Destroy()
   dev.off()

}

#############################################################################################
## in Exploring missing patterns menu
## missing rates per chip

#############################################################################################
### in Impute menu
## kNNimpute

Mnu.kNNimpute.handler<-function(w,u=NULL)
{   contents <- "Computing kNN imputation...\n"
     Result.Save<<-c(Result.Save,contents)
     gtkTextBufferInsertAtCursor(new.text,contents)
     gtkTextViewSetBuffer(outputText,new.text)
     subwindow("kNN")
     contents <- "...kNN imputation is DONE!\n"
     Result.Save<<-c(Result.Save,contents)
     gtkTextBufferInsertAtCursor(new.text,contents)
     gtkTextViewSetBuffer(outputText,new.text)


}
# LLS impute

Mnu.LLSimpute.handler<-function(w,u=NULL)
{   pc.n.flag<<-FALSE
    contents <- "Computing LLS imputation...\n"
    Result.Save<<-c(Result.Save,contents)
    gtkTextBufferInsertAtCursor(new.text,contents)
    gtkTextViewSetBuffer(outputText,new.text)
    subwindow("LLS")


}
## RLSP

Mnu.RLSPimpute.handler<-function(w,u=NULL)
{   pc.n.flag<<-FALSE
     contents <- "Computing RLSP imputation...\n"
     Result.Save<<-c(Result.Save,contents)
     gtkTextBufferInsertAtCursor(new.text,contents)
     gtkTextViewSetBuffer(outputText,new.text)
     subwindow("RLSP")


}

## BPCAimpute
Mnu.BPCAimpute.handler<-function(w,u=NULL)
{    contents <- "Computing BPCA imputation...\n Iteration will be up to 200\n"
     Result.Save<<-c(Result.Save,contents)
     gtkTextBufferInsertAtCursor(new.text,contents)
     gtkTextViewSetBuffer(outputText,new.text)
     subwindow("BPCA")
     contents <- "...BPCA imputation is DONE!\n"
     Result.Save<<-c(Result.Save,contents)
     gtkTextBufferInsertAtCursor(new.text,contents)
     gtkTextViewSetBuffer(outputText,new.text)


}

## User defined impute 1
Mnu.UserDef1.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open True Value file")
    ok <- win[["OkButton"]]
    cancel <- win[["CancelButton"]]
    ok$AddCallback("clicked", function(x) {
         User1.impute<<-read.csv(win$GetFilename(),header=T)
         rownames(User1.impute)<<-as.character(User1.impute[,1])
         User1.impute<<-User1.impute[,-1]
         if(sum(as.numeric( dim(User1.impute)==dim(data.missing)))!=2)
         { warning.gtk("imputed data is different from missing data")
         } else
        {  User1.impute<<-round(User1.impute[is.na(data.missing)],4) }
         missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),          kNN.impute,LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
         removeMissingList(nrow(missing.Name))
         showMissingList(missing.Name)
         win$Hide()
         win$Destroy()
    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })


}

## User defined impute 1
Mnu.UserDef2.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open True Value file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         User2.impute<<-read.csv(win$GetFilename(),header=T)
         rownames(User2.impute)<<-as.character(User2.impute[,1])
         User2.impute<<-User2.impute[,-1]
         if(sum(as.numeric( dim(User2.impute)==dim(data.missing)))!=2)
         { warning.gtk("imputed data is different from missing data")
         } else
        {  User2.impute<<-round(User2.impute[is.na(data.missing)],4) }
         missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute, 
                                         LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
         removeMissingList(nrow(missing.Name))
         showMissingList(missing.Name)
         win$Hide()
         win$Destroy()
    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })


}

subwindow<-function(method)
{  
   Method<<-method
   ok.handler<-function(w,u=NULL)
   {    K<<-as.numeric(gtkEntryGetText(KK))
        if(method=="kNN")
        { kNNimpute(K)
        }else if(method=="BPCA")
        { BPCAfill(data.missing)
        }else if(method=="LLS")
        { LSimpute(data.missing)
        }else if(method=="RLSP")
        { RLSPimpute(data.missing)
        }
        k.input$Hide()
        k.input$Destroy()
   }
  cancel.handler<-function(w,u=NULL)
   {
        k.input$Hide()
        k.input$Destroy()
   }
  check.handler<-function(w,u=NULL)
   {   win<- gtkFileSelection(paste("Open",Method, "imputed value file"))
       ok <- win[["OkButton"]]
       cancel <- win[["CancelButton"]]
       ok$AddCallback("clicked", function(x) {
       impute<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
       rownames(impute)<<-as.character(impute[,1])
       impute<<-impute[,-1]
       if(sum(as.numeric( dim(impute)==dim(data.missing)))!=2)
       {   warning.gtk("imputed data is different from missing data")
       } else
       {   missing<-impute[is.na(data.missing)]
           if(method=="kNN")
           {   kNN.impute<<-round(missing,4)
           }else if(method=="BPCA")
           {   BPCA.impute<<-round(missing,4)
           }else if(method=="LLS")
           {   LLS.impute<<-round(missing,4); 
           }else if(method=="RLSP")
           {   RLSP.impute<<-round(missing,4)
           }
           missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute,
                                           LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
           removeMissingList(nrow(missing.Name))
           showMissingList(missing.Name)
       }
       win$Hide()
       win$Destroy() 
    }) 
    cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
    })
    k.input$Hide()
    k.input$Destroy()
   }
   k.input <<- gtkWindow (show = FALSE) 
   Hpan1<-gtkVPanedNew(FALSE)
   Hpan2<-gtkVPanedNew(FALSE)
   box1 <- gtkHBox(FALSE, 3)
   box1$PackStart(gtkLabel("   K =  "))
   box1$PackEnd(KK<<-gtkEntry())
   gtkEntrySetText(KK,10)
   box1Frame <- gtkFrameNew("input K")
   gtkContainerAdd(box1Frame,box1)
   gtkWidgetShow(box1Frame)
   checkbox1<-gtkButton("read imputed value from file")
   checkbox1$SetUsize(200,25)
   gtkAddCallback(checkbox1, "clicked", check.handler)
   checkbox1Frame <- gtkFrameNew("input file option") 
   gtkContainerAdd(checkbox1Frame,checkbox1)
   gtkWidgetShow(checkbox1Frame)

   button.ok<-gtkButton("Ok")
   button.ok$SetUsize(50,25)
   button.cancel<-gtkButton("Cancel")
   button.cancel$SetUsize(50,25)
   gtkAddCallback(button.ok, "clicked", ok.handler)
   gtkAddCallback(button.cancel, "clicked", cancel.handler)
   Button<-gtkHBoxNew()
   gtkBoxPackStart(Button,button.ok,TRUE,FALSE,0)
   gtkBoxPackStart(Button,button.cancel,TRUE,FALSE,0)
   gtkPanedPack1(Hpan1,box1Frame,TRUE,TRUE)
   gtkPanedPack2(Hpan1,Button,TRUE,TRUE)
   Hpan1$Show()
   gtkPanedPack1(Hpan2,Hpan1,TRUE,TRUE)
   gtkPanedPack2(Hpan2,checkbox1Frame,TRUE,TRUE)
   Hpan2$Show()
   variable.box<-gtkVBoxNew(FALSE,0)
   gtkBoxPackStart(variable.box,Hpan2,TRUE,TRUE,0)
   k.input$Add(variable.box)
   k.input$Show()


}

#############################################
## imputation calculation
#
kNNimpute<-function(K)
{  kNNimpute.temp<<-impute.knn(as.matrix(data.missing),K)
   kNN.impute<<-round(kNNimpute.temp[is.na(data.missing)],4)
   missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute,
                                   LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
   removeMissingList(nrow(missing.Name))
   showMissingList(missing.Name)


}

LSimpute<-function(data.missing,mink=10)
{ messages<-"  imputing missing values \n using LS method \n\n  Do Not Close This window!\n
                This window will be closed automatically "
    warn.win <- gtkWindow(show = FALSE)
    box <- gtkVBox(TRUE, 0)
    box$SetUsize(-1, 80)
    box$PackStart(gtkLabel(messages))
    box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
    warn.win$Add(box)
    warn.win$Show()
    LLS.impute<<-LSimputeMain(as.matrix(data.missing))
    LLS.impute<<-round(LLS.impute,4)
    missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute,
                                  LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
    removeMissingList(nrow(missing.Name))
    showMissingList(missing.Name)
    LS.flag<<-TRUE
    warn.win$hide()
    warn.win$destroy()


}

RLSPimpute<-function(mink=10)
{   messages<-"  imputing missing values \n using RLSP method \n\n  Do Not Close This window!\n
                        This window will be closed automatically "
    warn.win <- gtkWindow(show = FALSE)
    box <- gtkVBox(TRUE, 0)
    box$SetUsize(-1, 80)
    box$PackStart(gtkLabel(messages))
    box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
    warn.win$Add(box)
    warn.win$Show()

    RLSP.impute<<-RLSPimputeMain(as.matrix(data.missing))
    RLSP.impute<<-round(RLSP.impute,4)

    missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute,
                                    LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
    removeMissingList(nrow(missing.Name))
    showMissingList(missing.Name)
    RLSP.flag<<-TRUE
    warn.win$hide()
    warn.win$destroy()

}

similargene <- function(E,i, missidxj, nomissidxj,gene.include)
{    mm1 <- 1
     mm2 <-length(gene.include)
     AA1 <- as.matrix(E[i, nomissidxj])
     BB1 <- as.matrix(E[gene.include, nomissidxj])
     E.temp<-E[gene.include,]
     Distance<-apply((t(AA1%*%matrix(rep(1,nrow(BB1)),nrow=1))-BB1)**2,1,sum)

     sorted <- sort.list(Distance,decreasing=TRUE)
     result <- list()      
     result$A <- E.temp[sorted,nomissidxj]
     result$B <- as.matrix(E.temp[sorted,missidxj])
     result$w <- as.matrix(E[i,nomissidxj] )
     result   
}        

warning.gtk <- function(messages) 
{ warn.win <- gtkWindow(show = FALSE)
    box <- gtkVBox(TRUE, 0)
    box$SetUsize(-1, 80)
    box$PackStart(gtkLabel(messages))
    box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
    box.small$PackStart(button <- gtkButton("Ok"), expand = TRUE, fill = FALSE)
    button$SetUsize(60, 25)
    button$AddCallback("clicked", function(x) {
      warn.win$Hide()
      warn.win$Destroy()
    })
    warn.win$Add(box)
    warn.win$Show()
 }

BPCAfill<-function(x999,k=NULL,maxepoch=NULL)
{  messages<-"  imputing missing values \n using BPCA method \n\n  Do Not Close This window!
                         \nThis window will be closed automatically  "
     warn.win <- gtkWindow(show = FALSE)
     box <- gtkVBox(TRUE, 0)
     box$SetUsize(-1, 80)
     box$PackStart(gtkLabel(messages))
     box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
     warn.win$Add(box)
     warn.win$Show()
     BPCA.impute<<-round(BPCAimputeMain(as.matrix(x999)),4)
     missing.Name<<-cbind(missing.Name.temp,round(AVE.impute,3),kNN.impute,
                                      LLS.impute,RLSP.impute,BPCA.impute,User1.impute,User2.impute)
     removeMissingList(nrow(missing.Name))
     showMissingList(missing.Name)
     warn.win$hide()
     warn.win$destroy()
 }


####
BPCA.dostep<-function(M,y)
{   q <- M$q;
    N <- M$N;
    d <- M$d;
    Rx <- diag(1,q)+M$tau*t(M$W)%*%M$W+M$SigW;
    Rxinv <- solve( Rx );
    idx <- c(M$gnomiss);
    n <- length(idx);
    dy <- y[idx,] - matrix(rep(1,n),nrow=n)%*%M$mu;
    x <- M$tau * Rxinv %*% t(M$W) %*% t(dy);
    T <- t(dy) %*% t(x);
    trS <- sum(sum(unlist(dy*dy)))
    for(n in 1:length(M$gmiss))
    {  i <- M$gmiss[n];
       dyo <- as.matrix(y[i,M$nomissidx[[i]]] - M$mu[M$nomissidx[[i]]])
       Wm <- matrix(M$W[M$missidx[[i]],],nrow=length(M$missidx[[i]]))
       Wo <- matrix(M$W[M$nomissidx[[i]],],nrow=length(M$nomissidx[[i]]))
       Rxinv <- solve( Rx - M$tau*t(Wm)%*%Wm );
       ex <- M$tau * t(Wo) %*% dyo;
       x <- Rxinv %*% ex;
       dym <- Wm %*% x;
       dy <- as.matrix(y[i,])
       dy[M$nomissidx[[i]]] <- dyo;
       dy[M$missidx[[i]]] <- dym;
       M$yest[i,] <- t(dy) + M$mu;
       T <- T + dy%*%t(x);
       T[M$missidx[[i]],] <- T[M$missidx[[i]],] + Wm %*% Rxinv;
       trS = trS + t(dy)%*%dy + 
            length(M$missidx[[i]])/M$tau +  sum(diag(Wm %*% Rxinv %*% t(Wm))) ;
    } 
   T = T/N;
   trS = trS/N;
   Rxinv <- solve(Rx);
   Dw <- Rxinv + M$tau*t(T)%*%M$W%*%Rxinv + diag(M$alpha)/N;
   Dwinv <- solve(Dw);
   M$W = T %*% Dwinv; 
   M$tau <- c(  (d+2*M$gtau0/N) / (trS-sum(diag((t(T)%*%M$W))) +
              (M$mu%*%t(M$mu)%*%M$gmu0+2*M$gtau0/M$btau0)/N))
   M$SigW <- Dwinv*(d/N);
   M$alpha <- (2*M$galpha0 + d)/
       	  (M$tau*diag(t(M$W)%*%M$W)+diag(M$SigW)+2*M$galpha0/M$balpha0);
  return(M)
}

###########
BPCA.initmodel<-function(y,q)
{ N<-nrow(y)
  d<-ncol(y)
  M<-list()
  M$N <- N;
  M$q <- q;
  M$d <- d;

  M$yest = y;
  M$missidx = list()
  M$nomissidx =  list()
  M$gnomiss = NULL;
  M$gmiss = NULL;
  for(i in 1:N)
  {
    missidx<- which(is.na(y[i,]));
    nomissidx = which(!is.na(y[i,]));
    M$missidx[[i]]<-missidx
    M$nomissidx[[i]]<-nomissidx
    if(length(missidx) == 0)
    {   M$gnomiss = cbind(M$gnomiss, i);
    } else
    {   M$gmiss = cbind(M$gmiss, i);
        M$yest[i,missidx] = 0;
    }
  }
  ynomiss <- y[M$gnomiss,];
  covy <- cov(M$yest);
  svd.result<-svd(covy)
  U<-svd.result$u;U<-U[,-c(d:d-(d-q-1))]
  S<-svd.result$d;S<-S[-c(d:d-(d-q-1))]
  V<-svd.result$v
  M$mu = matrix(0,nrow=1,ncol=d);
  for(j in 1:d)
  {  idx <- which(!is.na(y[,j]));
     M$mu[j] <- mean(y[idx,j]);
  }
  M$W  <- U %*% sqrt(diag(S));
  M$tau <- 1/( sum(diag(covy)) -sum(S) );
  taumax <- 1e10;
  taumin <- 1e-10;
  M$tau <- max( min( M$tau, taumax), taumin );
  M$galpha0 <- 1e-10;
  M$balpha0 <- 1;
  M$alpha <- (2*M$galpha0 + d)/(M$tau*diag(t(M$W)%*%M$W)+2*M$galpha0/M$balpha0);
  M$gmu0  <- 0.001;
  M$btau0 <- 1;
  M$gtau0 <- 1e-10;
  M$SigW <- diag(1,q);
  return(M)
}  


#############################################################################################
## in Comparison menu
## overall comparison
# t-test

Mnu.overall.compare.ttest.handler<-function(w,u=NULL)
{

 }


## MDS plot
Mnu.overall.compare.plot.handler<-function(w,u=NULL)
{ 

}


## Validation : compare to the true values - NRMSE
Mnu.comp.TrueValue.handler<-function(w,u=NULL)
{ 

}

## Validation : compare to the true values - NRMSE


Mnu.voting.handler<-function(w,u=NULL)
{ 

}

## anova
Mnu.ANOVA.handler<-function(w,u=NULL)
{

}

## imputation value comparison for each gene

Mnu.compare.plot.handler<-function(w,u=NULL)
{

}

## Guide to decide imputation method

Mnu.choose.impute.handler<-function(w,u=NULL)
{  
   

}

BPCAimputeMain<-function(x999)
{
 

}

RLSPimputeMain<-function(data.missing.temp,mink=10)
{
   

}



LSimputeMain<-function(data.missing.temp,mink=10)
{
   

}
############################################
############################################
# Main GUI
############################################
###########################################


 #-----  Main Window  -----#
   main <<- gtkWindow (show = FALSE) 
   main$SetTitle("arrayImpute")
   gtkWindowSetDefaultSize(main,970,450)
#   gtkWindowSetPolicy(main,TRUE,TRUE,FALSE)  #allow top level window resizable
   gtkWindowSetResizable(main,TRUE)
   main$SetUposition(40,100) #gtkWidgetSetUposition


# buttonNames<-c("kNNimpute"","BPCAimpute","LSimpute","LLSimpute","RLSPimpute","SVRimpute")

 #-----  Gene List  -----#     (Temporary)
   geneList <<- gtkCListNew(9)
   gtkCListSetColumnWidth(geneList,0,75)
   gtkCListSetColumnWidth(geneList,1,75)
   gtkCListSetColumnWidth(geneList,2,75)
   gtkCListSetColumnWidth(geneList,3,75)
   gtkCListSetColumnWidth(geneList,4,75)
   gtkCListSetColumnWidth(geneList,5,75)
   gtkCListSetColumnWidth(geneList,6,75)
   gtkCListSetColumnWidth(geneList,7,75)
   gtkCListSetColumnWidth(geneList,8,75)

   gtkCListColumnTitlesShow(geneList)
   label0 <- gtkLabelNew("Chip name")	
   gtkCListSetColumnWidget(geneList,0,label0)
   label1 <- gtkLabelNew("Gene name")	
   gtkCListSetColumnWidget(geneList,1,label1)
   label2 <- gtkLabelNew("average")	
   gtkCListSetColumnWidget(geneList,2,label2)
   label3 <- gtkLabelNew("kNNimpute")	
   gtkCListSetColumnWidget(geneList,3,label3)
   label4 <- gtkLabelNew("LLSimpute")	
   gtkCListSetColumnWidget(geneList,4,label4)
   label5 <- gtkLabelNew("RLSPimpute")	
   gtkCListSetColumnWidget(geneList,5,label5)
   label6 <- gtkLabelNew("BPCAimpute")	
   gtkCListSetColumnWidget(geneList,6,label6)
   label7 <- gtkLabelNew("UserDef.1")	
   gtkCListSetColumnWidget(geneList,7,label7)
   label8 <- gtkLabelNew("UserDef.2")	
   gtkCListSetColumnWidget(geneList,8,label8)
 
   gtkCListColumnTitlesPassive(geneList)
   gtkCListSetColumnResizeable(geneList,0,FALSE)
   gtkCListSetSelectionMode(geneList,GtkSelectionMode[4])


 #-----  Chip List  -----#     (Temporary)
   ChipList <<- gtkCListNew(1)
#   gtkCListSetColumnWidth(ChipList,0,85)

   gtkCListColumnTitlesShow(ChipList)
   label11 <- gtkLabelNew("Chip name")	
   gtkCListSetColumnWidget(ChipList,0,label11)
   gtkCListColumnTitlesPassive(ChipList)
   gtkCListSetColumnResizeable(ChipList,0,FALSE)
   gtkCListSetSelectionMode(ChipList,GtkSelectionMode[4])

 #----- output screen -----#
    outputText <<- gtkTextViewNew()
    gtkTextViewSetEditable(outputText,FALSE)
#    outputText$SetEditable(FALSE)
    new.text<<-gtkTextBufferNew()  

#----- Scroll Window 1&2 -----#
 
   scrollwindow1 <- gtkScrolledWindowNew()	# to go into the geneList frame
   scrollwindow1$SetUsize(700,70)
   gtkScrolledWindowSetPolicy(scrollwindow1,GtkPolicyType[2],GtkPolicyType[2])
   scrollwindow1$Add(geneList)
   gtkContainerSetBorderWidth(scrollwindow1,7)

   scrollwindow2 <- gtkScrolledWindowNew()	# to go into the geneList frame
   scrollwindow2$SetUsize(170,35)
   gtkScrolledWindowSetPolicy(scrollwindow2,GtkPolicyType[2],GtkPolicyType[2])
   scrollwindow2$Add(outputText)
   gtkContainerSetBorderWidth(scrollwindow2,3)

   scrollwindow3 <- gtkScrolledWindowNew()	# to go into the geneList frame
   scrollwindow3$SetUsize(100,25)
   gtkScrolledWindowSetPolicy(scrollwindow3,GtkPolicyType[2],GtkPolicyType[2])
   scrollwindow3$Add(ChipList)
   gtkContainerSetBorderWidth(scrollwindow3,3)

#-----------------------------------#
#               menu
#-----------------------------------#

   menu1 <- gtkMenuNew()
   menu2 <- gtkMenuNew()
   menu3 <- gtkMenuNew()
   menu4 <- gtkMenuNew()
   menu21 <- gtkMenuNew()
   menu41 <- gtkMenuNew()
   menu42 <- gtkMenuNew()
 
#----------

   show1 <- gtkMenuItemNewWithLabel("File")
   item11 <- gtkMenuItemNewWithLabel("Read Microarray data with missing")
   item12 <- gtkMenuItemNewWithLabel("Read Gene Location with pin information")
   item13 <- gtkMenuItemNewWithLabel("Save Result")
   item14 <- gtkMenuItemNewWithLabel("Exit")

   show3 <- gtkMenuItemNewWithLabel("Impute")
   item31 <-gtkMenuItemNewWithLabel("kNNimpute")
   item32 <-gtkMenuItemNewWithLabel("LLSimpute")
   item33 <-gtkMenuItemNewWithLabel("RLSPimpute")
   item34 <-gtkMenuItemNewWithLabel("BPCAimpute")
   item35 <-gtkMenuItemNewWithLabel("The User Defined impute 1")
   item36 <-gtkMenuItemNewWithLabel("The User Defined impute 2")
 
   show4 <- gtkMenuItemNewWithLabel("Comparison")   
   item41 <- gtkMenuItemNewWithLabel("Overall comparison")
   item411 <- gtkMenuItemNewWithLabel("t-test")
   item412 <- gtkMenuItemNewWithLabel("MDS plot")
   item42 <- gtkMenuItemNewWithLabel("Validation : Compare with True Values")
   item421 <- gtkMenuItemNewWithLabel("NRMSE")
   item422 <- gtkMenuItemNewWithLabel("Voting on imputation methods")
   item423 <- gtkMenuItemNewWithLabel("ANOVA")
   item43 <- gtkMenuItemNewWithLabel("Imputation Value Comparison for each gene")
   item44 <- gtkMenuItemNewWithLabel("Guide to decide imputation method")

#-----------

   gtkMenuShellAppend(menu1,item11)
   gtkMenuShellAppend(menu1,item12)
   gtkMenuShellAppend(menu1,item13)
   gtkMenuShellAppend(menu1,item14)
 
   gtkMenuShellAppend(menu3,item31) # setting parameters
   gtkMenuShellAppend(menu3,item32) # kNNimpute
   gtkMenuShellAppend(menu3,item33) # LLS/RLSP impute
   gtkMenuShellAppend(menu3,item34) # BPCA impute
   gtkMenuShellAppend(menu3,item35) # User defined impute
   gtkMenuShellAppend(menu3,item36) # User defined impute

   gtkMenuShellAppend(menu4,item41) #
   gtkMenuShellAppend(menu4,item42)
   gtkMenuShellAppend(menu4,item43)
   gtkMenuShellAppend(menu4,item44)

   gtkMenuShellAppend(menu41,item411) #
   gtkMenuShellAppend(menu41,item412) #
   gtkMenuShellAppend(menu42,item421)
   gtkMenuShellAppend(menu42,item422)
   gtkMenuShellAppend(menu42,item423)



#------------

   gtkWidgetShow(show1)
   gtkWidgetShow(item11) 
   gtkWidgetShow(item12) 
   gtkWidgetShow(item13) 
   gtkWidgetShow(item14) 
 
    gtkWidgetShow(show3)    
   gtkWidgetShow(item31) 
   gtkWidgetShow(item32)
   gtkWidgetShow(item33)
   gtkWidgetShow(item34)
   gtkWidgetShow(item35)   
   gtkWidgetShow(item36)   

   gtkWidgetShow(show4)    
   gtkWidgetShow(item41) 
   gtkWidgetShow(item411) 
   gtkWidgetShow(item412) 
   gtkWidgetShow(item42) 
   gtkWidgetShow(item43) 
   gtkWidgetShow(item44) 
   gtkWidgetShow(item421) 
   gtkWidgetShow(item422) 
   gtkWidgetShow(item423) 
  
#-----------

   gtkMenuItemSetSubmenu(show1,menu1)
   gtkMenuItemSetSubmenu(show3,menu3)  
   gtkMenuItemSetSubmenu(show4,menu4)  
   gtkMenuItemSetSubmenu(item41,menu41)
   gtkMenuItemSetSubmenu(item42,menu42)
#--------

   menuBar<-gtkMenuBarNew()
   gtkWidgetShow(menuBar)

   gtkMenuShellAppend(menuBar,show1)
   gtkMenuShellAppend(menuBar,show3)
   gtkMenuShellAppend(menuBar,show4)


#------ callback function

   gtkAddCallback(item11, "activate", Mnu.read.file.handler)
   gtkAddCallback(item12, "activate", Mnu.read.location.file.handler)
   gtkAddCallback(item13, "activate", Mnu.save.result.handler)
   gtkAddCallback(item14, "activate", Mnu.exit.handler)

   gtkAddCallback(item31, "activate", Mnu.kNNimpute.handler)
   gtkAddCallback(item32, "activate", Mnu.LLSimpute.handler)
   gtkAddCallback(item33, "activate", Mnu.RLSPimpute.handler)
   gtkAddCallback(item34, "activate", Mnu.BPCAimpute.handler)
   gtkAddCallback(item35, "activate", Mnu.UserDef1.handler)
   gtkAddCallback(item36, "activate", Mnu.UserDef2.handler)
 
   gtkAddCallback(item411, "activate", Mnu.overall.compare.ttest.handler)
   gtkAddCallback(item412, "activate", Mnu.overall.compare.plot.handler)
   gtkAddCallback(item421, "activate", Mnu.comp.TrueValue.handler)
   gtkAddCallback(item422, "activate", Mnu.voting.handler)
   gtkAddCallback(item423, "activate", Mnu.ANOVA.handler)
   gtkAddCallback(item43, "activate", Mnu.compare.plot.handler)
   gtkAddCallback(item44, "activate", Mnu.choose.impute.handler)
 
   gtkAddCallback(geneList, "select-row", geneListSelectRow.handler)
   gtkAddCallback(geneList, "unselect-row", geneListUnselectRow.handler)
   gtkAddCallback(ChipList, "select-row", chipListSelectRow.handler)
   gtkAddCallback(ChipList, "unselect-row", chipListUnselectRow.handler)

#------------------------------------------------------------------#
#------------------------------------------------------------------#


#-----  Frames  -----#  
#-----  Frame : GeneList List -----#
   GeneListFrame <- gtkFrameNew("Missing Imputation")
   gtkContainerAdd(GeneListFrame,scrollwindow1)
   gtkWidgetShow(GeneListFrame)


#----- Frame : OutputTable ---------#
   outputTableFrame<-gtkFrameNew("Notes")
   gtkContainerAdd(outputTableFrame,scrollwindow2)
   gtkWidgetShow(outputTableFrame)


#----- Frame : ChipList ---------#
   ChipListFrame<-gtkFrameNew("Chip List")
   gtkContainerAdd(ChipListFrame,scrollwindow3)
   gtkWidgetShow(ChipListFrame)


#----- Packing and Show -----#

    Hpan <-gtkHPanedNew(FALSE) # to put GeneListFrame & ChipListFrame together
    Vpan <-gtkVPanedNew(FALSE)
   gtkPanedPack1(Vpan,outputTableFrame,TRUE,TRUE) 
   gtkPanedPack2(Vpan,ChipListFrame,TRUE,TRUE) 
   Vpan$Show() 
   gtkPanedPack1(Hpan,Vpan,TRUE,TRUE)
   gtkPanedPack2(Hpan,GeneListFrame,TRUE,TRUE)
   Hpan$Show()
   mainLayout <- gtkVBoxNew(FALSE,0)
   gtkBoxPackStart(mainLayout,menuBar,FALSE,FALSE,0)
   gtkBoxPackStart(mainLayout,Hpan,TRUE,TRUE,0)
   main$Add(mainLayout)
   main$Show()

}
