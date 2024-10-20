# hg38vcf.to.chr.info.ok.awk

BEGIN{
    OFS="\t";
}

$1~/#/&&$1~/contig/{
    gsub(/chr/,"",$0); print;
}

$1~/#/&&$1!~/contig/{
    print;
}

$1!~/#/{
    gsub(/chr/,"",$1);
    $8=".";
    print;
}
