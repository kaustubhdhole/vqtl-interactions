
import csv
import numpy as np
import time
import sys

#######################################################
#             vQTL - vQTL Mapping             	      #
#######################################################
# Kaustubh Dhole || BITS Pilani, Goa || TIFR, Mumbai  #
#######################################################

fo = open('phenogeno.csv', 'r')
data = [[word.split(',') for word in line.split()][0] for line in fo]
rows = zip(*data)
fo.close()


pheno_data = {}
for row in rows[0:34]:
	pheno_data[row[0]]= np.array(row[3:],dtype='float')


geno_data = {}
for row in rows[34:]:
	geno_data[row[0]]= np.array(row[3:])

env=['Paraquat','Congo_red','Zeocin','Neomycin','Diamide','Magnesium_Sulfate','Tunicamycin','Indoleacetic_Acid','Formamide','Hydrogen_Peroxide','Hydroquinone','Hydroxyurea','x4NQO','Cadmium_Chloride','Cobalt_Chloride','Copper','Lithium_Chloride','Maltose','Calcium_Chloride','x4-Hydroxybenzaldehyde','SDS','Galactose','Mannose','x6-Azauracil','Magnesium_Chloride','Ethanol','Raffinose','Lactate','Lactose','Xylose','Sorbitol','Trehalose','x5-Fluorocytosine','x5-Fluorouracil']
env=['Paraquat','Congo_red','Zeocin','Neomycin','Diamide','Magnesium_Sulfate','Tunicamycin','Indoleacetic_Acid','Formamide','Hydrogen_Peroxide','Hydroquinone','Hydroxyurea','x4NQO','Cadmium_Chloride','Cobalt_Chloride','Copper','Lithium_Chloride','Maltose','Calcium_Chloride','x4.Hydroxybenzaldehyde','SDS','Galactose','Mannose','x6.Azauracil','Magnesium_Chloride','Ethanol','Raffinose','Lactate','Lactose','Xylose','Sorbitol','Trehalose','x5.Fluorocytosine','x5.Fluorouracil']

for e in env:
		# output file
		fw=open(r'Interaction'+'//'+'vQTL_'+e+'.csv','a+')
		fw.write('snp1'+','+'snp2'+','+'LODi'+','+'LODa'+','+'LODf'+','+'perm_value'+','+'mbb'+','+'mbr'+','+'mrr'+','+'mrb'+','+'vbb'+','+'vbr'+','+'vrr'+','+'vrb'+'\n')
		fw.close()

		fil=open(r'vQTL'+'//'+'e '+'.csv','r+');
		markers1 = [line.split(",")[1] for line in fil][2:]
		fil.close();
		for i in range(0,len(markers1)):
			markers1[i]=markers1[i].split("X")[1]

		markers1=[marker.split("\"")[0] for marker in markers1]

		markers = list(set(markers1))
		num_markers = len(markers); 
		num_comparisons = num_markers*(num_markers-1)/2.0
		print "no of comparisons=",num_comparisons
		condition = e
		#print "\nMarkers for condition ", condition, markers
		#print "\nNumber of markers: ", len(markers)

		print "Condition = ", condition
		pheno = np.copy(pheno_data[condition])
		numstrains = len(pheno)

		LOD_distribution = []
		num_iter = 10000
		for it in range(num_iter):
			count = 0
			LODmax = 0
			np.random.shuffle(pheno)
			#count = 0
			for i, snp1 in enumerate(markers):
				for j, snp2 in enumerate(markers):
					if i>j:
					
						g00 = np.logical_and(geno_data[snp1]=='BB',geno_data[snp2]=='BB')
						g10 = np.logical_and(geno_data[snp1]=='RR',geno_data[snp2]=='BB')
						g01 = np.logical_and(geno_data[snp1]=='BB',geno_data[snp2]=='RR')
						g11 = np.logical_and(geno_data[snp1]=='RR',geno_data[snp2]=='RR')

						#The following step is later used to exclude sites with missing genotypes at one of the two SNPs
						g_all = g00+g10+g01+g11

						#g_all = np.logical_or(g00,np.logical_or(g01,np.logical_or(g10,g11)))
						#print np.logical_and(g00,np.logical_and(g01,np.logical_and(g10,g11)))

						mdat = np.ma.masked_array(pheno,np.isnan(pheno))

						a = pheno[g00]
						mdata = np.ma.masked_array(a,np.isnan(a))
						b = pheno[g10]
						mdatb = np.ma.masked_array(b,np.isnan(b))
						c = pheno[g01]
						mdatc = np.ma.masked_array(c,np.isnan(c))
						d = pheno[g11]
						mdatd = np.ma.masked_array(d,np.isnan(d))

						n00 = np.size(a[~np.isnan(a)])
						n10 = np.size(b[~np.isnan(b)])
						n01 = np.size(c[~np.isnan(c)])
						n11 = np.size(d[~np.isnan(d)])

						n = (n00+n01+n10+n11)

						mdatA=abs(mdata-np.mean(mdata))
						mdatB=abs(mdatb-np.mean(mdatb))
						mdatC=abs(mdatc-np.mean(mdatc))
						mdatD=abs(mdatd-np.mean(mdatd))

						x00 = np.mean(mdatA)
						x10 = np.mean(mdatB)
						x01 = np.mean(mdatC)
						x11 = np.mean(mdatD)

						overall_mean = (n00*x00+n10*x10+n01*x01+n11*x11)/(n00+n01+n10+n11)
						#overall_mean = np.mean(mdat)
						#zpheno=pheno

						#zpheno[g00]=abs(pheno[g00]-np.mean(mdata))
						#zpheno[g01]=abs(pheno[g01]-np.mean(mdatb))
						#zpheno[g10]=abs(pheno[g10]-np.mean(mdatc))
						#zpheno[g11]=abs(pheno[g11]-np.mean(mdatd))
						#[chr1,pos1] = snp1.split('_')
						#[chr2,pos2] = snp2.split('_')
						chr1=snp1.split('_')[1]; pos1=snp1.split('_')[2]; 
						chr2=snp2.split('_')[1]; pos2=snp2.split('_')[2]; 


						if not n00*n01*n10*n11 == 0 and n > 0.9*numstrains and (not (chr1 == chr2) or ((chr1 == chr2) and np.abs(int(pos1)-int(pos2)) > 10000)): 
							count += 1

							#Log Likelihood of data fitting a single mean
							#a = np.square(pheno - overall_mean)[g_all]
							#a = np.square(zpheno - overall_mean)[g_all]
							a1 = np.square(mdatA - overall_mean);mdata1 = np.ma.masked_array(a1,np.isnan(a1));
							a2 = np.square(mdatB - overall_mean);mdata2 = np.ma.masked_array(a2,np.isnan(a2));
							a3 = np.square(mdatC - overall_mean);mdata3 = np.ma.masked_array(a3,np.isnan(a3));
							a4 = np.square(mdatD - overall_mean);mdata4 = np.ma.masked_array(a4,np.isnan(a4));
							
							RSS0 = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

							#Log Likelihood of data fitting 4 means
							mu = x00
							beta1 = x10 - mu
							beta2 = x01 - mu
							gamma = x11 - mu - beta1 - beta2
							#print mu, mu+beta1, mu+beta2, mu+beta1+beta2+gamma
				
							#a = np.square(zpheno - mu*g00 - (mu + beta1)*g10 - (mu + beta2)*g01 - (mu + beta1 + beta2 + gamma)*g11)[g_all]
							#mdat = np.ma.masked_array(a,np.isnan(a))
							a1 = np.square(mdatA - mu); mdata1 = np.ma.masked_array(a1,np.isnan(a1))
							a2 = np.square(mdatB - (mu + beta1)); mdata2 = np.ma.masked_array(a2,np.isnan(a2))
							a3 = np.square(mdatC - (mu + beta2)); mdata3 = np.ma.masked_array(a3,np.isnan(a3))
							a4 = np.square(mdatD - (mu + beta1 + beta2 + gamma)); mdata4 = np.ma.masked_array(a4,np.isnan(a4))
							RSSF = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

							#Log likelihood of data fitting additive model
							den = n01*n10*(n00+n11) + n00*n11*(n01 + n10)
							mu = (n00*(n10*n11 + n01*(n10 + n11))*x00 + n01*n10*n11*(x01 + x10 - x11))/den
							beta1 = ((n00 + n10)*n01*n11*(x11-x01) + n00*n10*(n01 + n11)*(x10-x00))/den
							beta2 = ((n00 + n01)*n10*n11*(x11-x10) + n00*n01*(n10 + n11)*(x01-x00))/den

							#beta1 = (np.mean(mdatb) + np.mean(mdatd) - np.mean(mdatc) - np.mean(mdata))/2.0
							#beta2 = (np.mean(mdatc) + np.mean(mdatd) - np.mean(mdatb) - np.mean(mdata))/2.0

							#np.square
							a1 = np.square(mdatA - mu); mdata1 = np.ma.masked_array(a1,np.isnan(a1));
							a2 = np.square(mdatB - (mu + beta1)); mdata2 = np.ma.masked_array(a2,np.isnan(a2));
							a3 = np.square(mdatC - (mu + beta2)); mdata3 = np.ma.masked_array(a3,np.isnan(a3));
							a4 = np.square(mdatD - (mu + beta1 + beta2)); mdata4 = np.ma.masked_array(a4,np.isnan(a4));
							
							RSSA = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

							LODi = -n/2.0 * np.log10(RSSF/RSSA)
							#print snp1, snp2, LODi
							#print LOD_list.append(LODi)
							if LODi > LODmax:
								LODmax = LODi
								SNPmax = [snp1,snp2]
			
			#print "max: ", LODmax, count
			print('Iteration :'+ str(it)+' '+e)	
			LOD_distribution.append(LODmax)
		#print sorted(LOD_distribution)
#lodfile=open('lodfile_34_paraquat.csv','w+')#____revert_____#
#for i in range(num_iter):
#	lodfile.write(str(LOD_distribution[i])+','+'\n')
#lodfile.close()


		#NOW FOR THE REAL DEAL
		pheno = np.copy(pheno_data[condition])
		numstrains = len(pheno)

		LOD = np.array([[0.0 for i in range(len(markers))] for i in range(len(markers))])
		
		#count = 0
		for i, snp1 in enumerate(markers):
			for j, snp2 in enumerate(markers):
				if i>j:
				#if snp1!=snp2:

			#snp1 = 'chr04_29152'
			#snp2 = 'chr10_404421'
					
					g00 = np.logical_and(geno_data[snp1]=='BB',geno_data[snp2]=='BB')
					g10 = np.logical_and(geno_data[snp1]=='RR',geno_data[snp2]=='BB')
					g01 = np.logical_and(geno_data[snp1]=='BB',geno_data[snp2]=='RR')
					g11 = np.logical_and(geno_data[snp1]=='RR',geno_data[snp2]=='RR')
					#The following step is later used to exclude sites with missing genotypes at one of the two SNPs
					g_all = g00+g10+g01+g11

					#g_all = np.logical_or(g00,np.logical_or(g01,np.logical_or(g10,g11)))
					#print np.logical_and(g00,np.logical_and(g01,np.logical_and(g10,g11)))

					mdat = np.ma.masked_array(pheno,np.isnan(pheno))

					a = pheno[g00]
					mdata = np.ma.masked_array(a,np.isnan(a))
					b = pheno[g10]
					mdatb = np.ma.masked_array(b,np.isnan(b))
					c = pheno[g01]
					mdatc = np.ma.masked_array(c,np.isnan(c))
					d = pheno[g11]
					mdatd = np.ma.masked_array(d,np.isnan(d))

					n00 = np.size(a[~np.isnan(a)])
					n10 = np.size(b[~np.isnan(b)])
					n01 = np.size(c[~np.isnan(c)])
					n11 = np.size(d[~np.isnan(d)])

					n = (n00+n01+n10+n11)

					mdatA=abs(mdata-np.mean(mdata))
					mdatB=abs(mdatb-np.mean(mdatb))
					mdatC=abs(mdatc-np.mean(mdatc))
					mdatD=abs(mdatd-np.mean(mdatd))
					
					#zpheno=pheno

					#zpheno[g00]=abs(pheno[g00]-np.mean(mdata))
					#zpheno[g01]=abs(pheno[g01]-np.mean(mdatb))
					#zpheno[g10]=abs(pheno[g10]-np.mean(mdatc))
					#zpheno[g11]=abs(pheno[g11]-np.mean(mdatd))


					#bb1=zpheno[geno_data[snp1]=='BB'];bb1=np.ma.masked_array(bb1,np.isnan(bb1));m1bb=np.mean(bb1);v1bb=np.var(bb1)
					#rr1=zpheno[geno_data[snp1]=='RR'];rr1=np.ma.masked_array(rr1,np.isnan(rr1));m1rr=np.mean(rr1);v1rr=np.var(rr1)
					#bb2=zpheno[geno_data[snp2]=='BB'];bb2=np.ma.masked_array(bb2,np.isnan(bb2));m2bb=np.mean(bb2);v2bb=np.var(bb2)
					#rr2=zpheno[geno_data[snp2]=='RR'];rr2=np.ma.masked_array(rr2,np.isnan(rr2));m2rr=np.mean(rr2);v2rr=np.var(rr2)
					
					x00 = np.mean(mdatA)
					x10 = np.mean(mdatB)
					x01 = np.mean(mdatC)
					x11 = np.mean(mdatD)

					y00 = np.var(mdatA)
					y10 = np.var(mdatB)
					y01 = np.var(mdatC)
					y11 = np.var(mdatD)

					overall_mean = (n00*x00+n10*x10+n01*x01+n11*x11)/(n00+n01+n10+n11)
					#overall_mean = np.mean(mdat)
					
					#[chr1,pos1] = snp1.split('_')
					#[chr2,pos2] = snp2.split('_')
					chr1=snp1.split('_')[1]; pos1=snp1.split('_')[2]; 
					chr2=snp2.split('_')[1]; pos2=snp2.split('_')[2]; 


					if not n00*n01*n10*n11 == 0 and n > 0.9*numstrains and (not (chr1 == chr2) or ((chr1 == chr2) and np.abs(int(pos1)-int(pos2)) > 10000)): 
						#count += 1

						#Log Likelihood of data fitting a single mean
						#a = np.square(pheno - overall_mean)[g_all]
						#a = np.square(zpheno - overall_mean)[g_all]
						a1 = np.square(mdatA - overall_mean);mdata1 = np.ma.masked_array(a1,np.isnan(a1));
						a2 = np.square(mdatB - overall_mean);mdata2 = np.ma.masked_array(a2,np.isnan(a2));
						a3 = np.square(mdatC - overall_mean);mdata3 = np.ma.masked_array(a3,np.isnan(a3));
						a4 = np.square(mdatD - overall_mean);mdata4 = np.ma.masked_array(a4,np.isnan(a4));
							
						RSS0 = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

						#Log Likelihood of data fitting 4 means
						mu = x00
						beta1 = x10 - mu
						beta2 = x01 - mu
						gamma = x11 - mu - beta1 - beta2
						#print mu, mu+beta1, mu+beta2, mu+beta1+beta2+gamma
										#a = np.square(zpheno - mu*g00 - (mu + beta1)*g10 - (mu + beta2)*g01 - (mu + beta1 + beta2 + gamma)*g11)[g_all]
						mdat = np.ma.masked_array(a,np.isnan(a))
						a1 = np.square(mdatA - mu); mdata1 = np.ma.masked_array(a1,np.isnan(a1))
						a2 = np.square(mdatB - (mu + beta1)); mdata2 = np.ma.masked_array(a2,np.isnan(a2))
						a3 = np.square(mdatC - (mu + beta2)); mdata3 = np.ma.masked_array(a3,np.isnan(a3))
						a4 = np.square(mdatD - (mu + beta1 + beta2 + gamma)); mdata4 = np.ma.masked_array(a4,np.isnan(a4))
						RSSF = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

						#Log likelihood of data fitting additive model
						den = n01*n10*(n00+n11) + n00*n11*(n01 + n10)
						mu = (n00*(n10*n11 + n01*(n10 + n11))*x00 + n01*n10*n11*(x01 + x10 - x11))/den
						beta1 = ((n00 + n10)*n01*n11*(x11-x01) + n00*n10*(n01 + n11)*(x10-x00))/den
						beta2 = ((n00 + n01)*n10*n11*(x11-x10) + n00*n01*(n10 + n11)*(x01-x00))/den

						#beta1 = (np.mean(mdatb) + np.mean(mdatd) - np.mean(mdatc) - np.mean(mdata))/2.0
						#beta2 = (np.mean(mdatc) + np.mean(mdatd) - np.mean(mdatb) - np.mean(mdata))/2.0

						#np.square
						a1 = np.square(mdatA - mu); mdata1 = np.ma.masked_array(a1,np.isnan(a1));
						a2 = np.square(mdatB - (mu + beta1)); mdata2 = np.ma.masked_array(a2,np.isnan(a2));
						a3 = np.square(mdatC - (mu + beta2)); mdata3 = np.ma.masked_array(a3,np.isnan(a3));
						a4 = np.square(mdatD - (mu + beta1 + beta2)); mdata4 = np.ma.masked_array(a4,np.isnan(a4));
								
						RSSA = np.sum(mdata1) + np.sum(mdata2) +np.sum(mdata3) +np.sum(mdata4)

						#n = np.size(a[~np.isnan(a)])

						LODf = -n/2.0 * np.log10(RSSF/RSS0)
						LODa = -n/2.0 * np.log10(RSSA/RSS0)
						LODi = -n/2.0 * np.log10(RSSF/RSSA)
						fstat = (RSSA/RSSF - 1)*(n-2)
						#varf = (1 - RSSF/RSS0)
						#vara = (1 - RSSA/RSS0)
						#vari = varf - vara
						#print "F statistic = ", round((10**(2*LOD/n) - 1)*(n-2),3)
						
						#pval = 1-f.cdf(fstat,2-1,n-2)
						
						#pval_corrected_N_choose_2 = 1-((1-pval)**num_comparisons)
						#pval_corrected_N = 1-((1-pval)**num_markers)

						
						pval_perm = np.sum(np.array(LOD_distribution)>LODi)/float(num_iter)

						#if LODi >= 3:
						#if pval < 0.0005:
						if pval_perm <= 1:
							print ""
							print snp1, snp2
							#print "sample size: ",n
							#print n00, n10, n01, n11
							#print x00, x10, x01, x11
							print "LOD_f\tLOD_a\tLOD_i\tF\tp-val(uncorrected)\tp-val(corrected by N)\tp-val(corrected by N choose 2)\tp-val(permutation)"
							#print str(round(LODf,3))+"\t"+str(round(LODa,3))+"\t"+str(round(LODi,3))+"\t"+str(round(fstat,3))+"\t"+str(pval)+"\t"+str(round(pval_corrected_N,3))+"\t\t\t"+str(round(pval_corrected_N_choose_2,3))+"\t\t\t\t"+str(round(pval_perm,3))
							print str(round(LODf,3))+"\t"+str(round(LODa,3))+"\t"+str(round(LODi,3))+"\t"+str(round(fstat,3))+"\t\t\t"+str(round(pval_perm,3))
							print "LOD_f = ", round(LODf,3)
							print "LOD_a = ", round(LODa,3)
							print "LOD_i = ", round(LODi,3)
							print "F statistic = ", round(fstat,3)
						#	print "p-value = ", pval
							print(type(snp1));
							print(type(snp2));
							print(type(LODi));
							print(type(LODa));
							print(type(LODf));
							fw=open(r'Interaction'+'//'+'vQTL_'+e+'.csv','a+')
							#fw.write(snp1+','+snp2+','+LODi+','+LODa+','+LODf+','+pval+'\n')
							fw.write(snp1+','+snp2+','+str(LODi)+','+str(LODa)+','+str(LODf)+','+str(pval_perm)+','+str(x00)+','+str(x01)+','+str(x11)+','+str(x10)+','+str(y00)+','+str(y01)+','+str(y11)+','+str(y10)+'\n')
							fw.close()

						LOD[i][j] = LODi
						LOD[j][i] = LODi


			#plt.clf()
			#plt.subplots_adjust(bottom=0)
			#plt.imshow(LOD,interpolation='nearest', cmap=plt.cm.jet,aspect='auto', origin='lower',vmin=0,vmax=4)
			#plt.colorbar()
			#plt.yticks(range(len(markers)), markers, fontsize=4)
			#plt.xticks(range(len(markers)), markers, rotation='vertical',fontsize=4)
			#plt.savefig("heatmap"+"YP"+str(e)+"_"+str(p)+".png",bbox_inches='tight',dpi=900)
                        

