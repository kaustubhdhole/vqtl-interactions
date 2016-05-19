# Choose Environments here.
print("Specify the environments in interaction.R");
envnames = c("Paraquat","Congo_red","Zeocin","Neomycin","Diamide","Magnesium_Sulfate","Tunicamycin","Indoleacetic_Acid","Formamide","Hydrogen_Peroxide","Hydroquinone","Hydroxyurea","x4NQO","Cadmium_Chloride","Cobalt_Chloride","Copper","Lithium_Chloride","Maltose","Calcium_Chloride","x4.Hydroxybenzaldehyde","SDS","Galactose","Mannose","x6.Azauracil","Magnesium_Chloride","Ethanol","Raffinose","Lactate","Lactose","Xylose","Sorbitol","Trehalose","x5.Fluorocytosine","x5.Fluorouracil")  

DataFile3 = NULL;
for(e in envnames)
{ qtlfile = read.csv(file.path("2qtl",paste("Output_with_pval_",e,".csv",sep="")),stringsAsFactors=FALSE);
 vqtlfile = read.csv(file.path("2vqtl",paste("Output_with_permval_",e,".csv",sep="")),stringsAsFactors=FALSE);
 if(nrow(qtlfile)>0)
 {for(r in 1:nrow(qtlfile))
	{ 
	DataFile3 = rbind(DataFile3,c(e,qtlfile$snp1[r],qtlfile$snp2[r],"QTL-QTL","-",qtlfile$LODi[r],qtlfile$LODa[r],qtlfile$LODf[r]))
	}
 }
 if(nrow(vqtlfile)){
for(r in 1:nrow(vqtlfile))
	{ 
	DataFile3 = rbind(DataFile3,c(e,vqtlfile$snp1[r],vqtlfile$snp2[r],"-","vQTL-vQTL",vqtlfile$LODi[r],vqtlfile$LODa[r],vqtlfile$LODf[r]))
	}
	}
}
colnames(DataFile3)=c("Environment","M1","M2","QTL-QTL Interaction","vQTL-vQTL Interaction","LODi","LODa","LODf")
write.csv(DataFile3,"DataFile3.csv")
