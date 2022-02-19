awk '{n+=$1;N+=$2;printf "%s %s %.2f%%\n",$1,$2,$2==0?0:$1*100/$2}END{if (NR>1) printf "%s %s %.2f%%\n",n,N,N==0?0:n*100/N}'
