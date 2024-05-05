library(rhdf5)
library(Matrix)

set.seed(504)

RNAseqREF <- H5Fopen("HumanFetalBrainPool.h5")

h5ls(RNAseqREF)


region <- RNAseqREF$"/shoji/Region"
table(region)
dim(region)


telencephalon <- which(region== "Telencephalon")
head(telencephalon)

chunk1<- telencephalon[1:76355]
head(chunk1)
tail(chunk1)
dim(chunk1)

chunk2<- telencephalon[76356:152710]

chunk3<- telencephalon[152711:229065]

chunk4<- telencephalon[229066:305420]

chunk5<- telencephalon[305421:381776]

RNAseq.telencephalon.1 <- h5read(file= RNAseqREF, name = "/shoji/Expression", index = list(1:59480,chunk1),
)

RNAseq.telencephalon.2 <- h5read(file= RNAseqREF, name = "/shoji/Expression", index = list(1:59480,chunk2),
)

RNAseq.telencephalon.3 <- h5read(file= RNAseqREF, name = "/shoji/Expression", index = list(1:59480,chunk3),
)

RNAseq.telencephalon.4 <- h5read(file= RNAseqREF, name = "/shoji/Expression", index = list(1:59480,chunk4),
)

RNAseq.telencephalon.5 <- h5read(file= RNAseqREF, name = "/shoji/Expression", index = list(1:59480,chunk5),
)
RNAseq.telencephalon.1<- as(RNAseq.telencephalon.1, "dgCMatrix")
RNAseq.telencephalon.2<- as(RNAseq.telencephalon.2, "dgCMatrix")
RNAseq.telencephalon.3<- as(RNAseq.telencephalon.3, "dgCMatrix")
RNAseq.telencephalon.4<- as(RNAseq.telencephalon.4, "dgCMatrix")
RNAseq.telencephalon.5<- as(RNAseq.telencephalon.5, "dgCMatrix")

RNAseq.telencephalon <- cbind(RNAseq.telencephalon.1, RNAseq.telencephalon.2, RNAseq.telencephalon.3, RNAseq.telencephalon.4, RNAseq.telencephalon.5)

dim(RNAseq.telencephalon)

colnames(RNAseq.telencephalon)<- RNAseqREF$"/shoji/CellID"[telencephalon]

rownames(RNAseq.telencephalon)<-RNAseqREF$"/shoji/Gene"

tail(RNAseqREF$"/shoji/CellID", 10)

Age<- RNAseqREF$"/shoji/Age" [telencephalon]

#loading CellClass data
CellClass<- RNAseqREF$"/shoji/CellClass"[telencephalon]

table(CellClass)
dim(CellClass)

#loading CellCycleFraction data
CellCycleFraction<- RNAseqREF$"/shoji/CellCycleFraction" [telencephalon]

#loading Chemistry data
Chemistry<- RNAseqREF$"/shoji/Chemistry" [telencephalon]

Clusters<- RNAseqREF$"/shoji/Clusters"[telencephalon]

#loading donor data
Donor<- RNAseqREF$"/shoji/Donor" [telencephalon]

#loading NGenes data
NGenes<- RNAseqREF$"/shoji/NGenes" [telencephalon]

#loading Region data
Region<- RNAseqREF$"/shoji/Region" [telencephalon]

#loading SampleID
SampleID<- RNAseqREF$"/shoji/SampleID" [telencephalon]

#loading sex data
Sex<-RNAseqREF$"/shoji/Sex" [telencephalon]

#loading Subdivision data
Subdivision<- RNAseqREF$"/shoji/Subdivision" [telencephalon] 

#loading Subregion data
Subregion<- RNAseqREF$"/shoji/Subregion" [telencephalon]

#loading Tissue data
Tissue<- RNAseqREF$"/shoji/Tissue" [telencephalon]

AnnotationName <- RNAseqREF$"/shoji/AnnotationName" [telencephalon]

Class <- RNAseqREF$"/shoji/Class" [telencephalon]

table(Class)

#loading TopLevelCluster data
TopLevelCluster<- RNAseqREF$"/shoji/TopLevelCluster" [telencephalon]

#loading TotalUMIs data
TotalUMIs<- RNAseqREF$"/shoji/TotalUMIs" [telencephalon]

CellID<- RNAseqREF$"/shoji/CellID" [telencephalon]

RNAseq.meta.data<- data.frame(Clusters=Clusters, AnnotationName= AnnotationName, Class=Class, CellClass=CellClass, Age=Age,CellCycleFraction=CellCycleFraction,Chemistry=Chemistry,Donor=Donor,NGenes=NGenes,Region=Region,SampleID=SampleID, Sex=Sex,Subdivision=Subdivision, Subregion=Subregion, Tissue=Tissue,TopLevelCluster=TopLevelCluster,TotalUMIs=TotalUMIs, row.names = CellID)

library(dplyr)
library(Seurat)

RNAseq.seurat <- CreateSeuratObject(counts = RNAseq.telencephalon, project = "first.trimester", assay="RNA", meta.data = RNAseq.meta.data)

RNAseq.seurat
head(RNAseq.seurat)

saveRDS(RNAseq.seurat, "first_trimemster.rds")

