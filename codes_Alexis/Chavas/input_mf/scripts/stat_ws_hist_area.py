import numpy as np
import matplotlib.pyplot as plt

z=6
z=str(z)
a=np.loadtxt('/net/merzhin1/vol/safo/morgane/data/stat_ws_result_all.txt')
b=np.arange(16)

bins=[x + 1 for x in range(0, 5)]
diff1=[]
diff2=[]
diff3=[]
std1=[]
std2=[]
std3=[]
nb1=[]
nb2=[]
nb3=[]
for i in b[0:15:3]:
	diff1=np.append(a[i,0],diff1)
	diff2=np.append(a[i+1,0],diff2)
	diff3=np.append(a[i+2,0],diff3)
	std1=np.append(a[i,1],std1) 
	std2=np.append(a[i+1,1],std2)
	std3=np.append(a[i+2,1],std3)
	nb1=np.append(a[i,3],nb1)
	nb2=np.append(a[i+1,3],nb2)
	nb3=np.append(a[i+2,3],nb3)

diff2=diff2[::-1]	
diff1=diff1[::-1]	
diff3=diff3[::-1]
std1=std1[::-1]	
std2=std2[::-1]	
std3=std3[::-1]	
nb1=nb1[::-1]
nb2=nb2[::-1]
nb3=nb3[::-1]

s=200


fig=figure(figsize(25,17))
fig.suptitle('Statistique par catégorie de cyclone' ,fontsize = 24)
ax1 = fig.add_axes((0.05,0.10,0.40,0.80))
ax2 = fig.add_axes((0.55,0.10,0.40,0.80))

for j in range(len(bins)):
	ax1.bar(bins[j]-0.25,diff1[j],align='center',width=0.25,color='orange',label='SMOS' if j == 0 else "")  
	ax3=ax1.twinx()
	ax3.scatter( bins[j]-0.25,nb1[j],s=s,marker='o',color='orange',edgecolors='k')
	ax3.set_ylim([0,12500])
	ax3.set_ylabel('Nombre de points',fontsize=18)
	ax3.tick_params(axis = 'both', labelsize = 16)
	ax1.bar(bins[j],diff2[j],align='center',width=0.25,color='green',label='SMAP(RSS)' if j == 0 else "")  	

	ax3.scatter(bins[j],nb2[j],s=s,marker='o',color='green',edgecolors='k')
	ax3.set_ylim([0,12500])
	ax1.bar(bins[j]+0.25,diff3[j],align='center',width=0.25,color='blue',label='SMAP(JPL)' if j == 0 else "")

	ax3.scatter( bins[j]+0.25,nb3[j],s=s,marker='o',color='blue',edgecolors='k')
	ax3.set_ylim([0,12500])

	ax2.bar(bins[j]-0.25,std1[j],align='center',width=0.25,color='orange',label='SMOS' if j == 0 else "")  
	ax2.bar(bins[j],std2[j],align='center',width=0.25,color='green',label='SMAP(RSS)' if j == 0 else "")  
	ax2.bar(bins[j]+0.25,std3[j],align='center',width=0.25,color='blue',label='SMAP(JPL)' if j == 0 else "")

ax1.set_xlabel('Catégorie',fontsize=18)
ax1.set_ylabel('Biais',fontsize=18)
ax1.set_title('Biais',fontsize=22)
ax3.set_ylabel('Nombre de points',fontsize=18)
ax1.legend(fontsize=16)
ax1.grid()
ax1.set_ylim([-3.5,1.5])
ax1.tick_params(axis = 'both', labelsize = 16)
ax2.set_xlabel('Catégorie',fontsize=18)
ax2.set_ylabel('Std',fontsize=18)
ax2.set_title('Ecart-type' ,fontsize=22)
ax2.legend(fontsize=16)
ax2.grid()
ax2.set_ylim([0,8])
ax2.tick_params(axis = 'both', labelsize = 16)
savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/STAT_WS/stat_ws_cat_all_area.png')
plt.close('all')
