BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        split($0,a,"\t"); 
        para[a[1]]=1; 
        para[a[2]]=1; 
        agecl[a[1]]=a[10]; 
        agecl[a[2]]=a[10];
    }
} 

{
    split($1,a,">"); 
    split($4,b,":"); 
    split($5,c,":"); 
    bt[b[2]]=c[2]; 
    nbtr[b[2]]++; 
    trlist[b[2]]=(trlist[b[2]])(a[2])(",");
} 
    
END{
    OFS="\t"; 
    for(g in bt)
    {
        split(g,a,"."); 
        print g, bt[g], (para[a[1]]==1 ? 1 : 0), (para[a[1]]==1 ? cl(agecl[a[1]]) : "NA"), trlist[g], nbtr[g],
    }
} 

function cl(x){
    return (x=="young" ? "1.young" : (x=="old" ? "2.old" : ("3.veryold")));
} 