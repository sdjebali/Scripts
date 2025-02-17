# zotero2hal.awk
################
# The idea is to automate the submission of papers from our team (specially clinical papers
# that are many and with a lot of authors) into hal for the hceres 2025

# This script starts from:
##########################
# - a basename for the xml files to be generated and named base
# - three input files:
#   * fileRef1 is a tsv file with forename and surname of the authors for which we want to add the irsd affiliation
#     (at least one author per paper needs to be affiliated)
#   * fileRef2 is a template xml file to submit a paper in hal, with 5 places that should be replaced with content
#     - {TITLE} at two places with a string each time (the same one)
#     - {AUTHORS} at one place with several blocks of authors
#     - {JOURNAL} at one place with a block
#     - {DOI} at one place with a string
#   * main input file that is the zotero output in xml format with several papers that we want to submit to hal
# This script will produce:
###########################
# - as main output file a txt file with all the papers written one after the other in this form
#   << - title, authors, journal, year, idno, volume, issue, page, doi >>
# - as many hal compatible (TEI) xml files as we have papers in the main input file named file.no.xml where no goes from 1 to nb articles
# Note that the xml files that are output have many spaces before the information but it does not seem to matter for hal
# since it does matter for the main txt output file, we need to do this to it to make it right
# sed 's/  */ /g' zotero.output.txt > zotero.output.ok.txt
# - Predicting survival in patients with ‘non-high-risk’ acute variceal bleeding receiving β-blockers+ligation to prevent re-bleeding, Balcar L, Mandorfer M, Hernández-Gea V, Procopet B, Meyer E, Giráldez Á, Amitrano L, Villanueva C, Thabut D, Samaniego L, Silva-Junior G, Martinez J, Genescà J, Bureau C, Trebicka J, Herrera E, Laleman W, Palazón Azorín J, Alonso J, Gluud L, Ferreira C, Cañete N, Rodríguez M, Ferlitsch A, Mundi J, Grønbæk H, Hernandez Guerra M, Sassatelli R, Dell'Era A, Senzolo M, Abraldes J, Romero-Gómez M, Zipprich A, Casas M, Masnou H, Primignani M, Krag A, Nevens F, Calleja J, Jansen C, Catalina M, Albillos A, Rudler M, Tapias E, Guardascione M, Tantau M, Schwarzer R, Reiberger T, Laursen S, Lopez-Gomez M, Cachero A, Ferrarese A, Ripoll C, La Mura V, Bosch J, García-Pagán J, Journal of Hepatology, 2024, 01688278, 80, 1, 73-81, 10.1016/j.jhep.2023.10.007
# - Incidence and factors predictive of recurrent thrombosis in people with non-cirrhotic portal vein thrombosis, Baiges A, Procopet B, Silva-Junior G, Llop E, Tellez L, Darnell A, Garcia-Criado Á, Turon F, Nicoara-Farcau O, González-Alayón C, Larrue H, Magaz M, Olivas P, Perez-Campuzano V, Calleja J, Albillos A, Reverter J, Bureau C, Bosch J, Hernández-Gea V, Garcia-Pagán J, Journal of Hepatology, 2023, 01688278, 78, 1, 114-122, 10.1016/j.jhep.2022.08.023
# 1 (66 fields)
# 1 (139 fields)


# The idea is then to submit all of them one by one using a command like this one
# curl -X POST -d @art1.works.xml -v -u ${HAL_USERNAME}:${HAL_PASSWORD} https://api.archives-ouvertes.fr/sword/hal/ -H "Packaging:http://purl.org/net/sword-types/AOfr" -H "Content-Type:text/xml"
# where art1.works.xml is one of the xml file produced by this script

# Several assumptions are done here:
####################################
# - the title of the article is in the row that follows the row <title level="a">
# - the doi of the article is in the row that follows the row <idno type="DOI">
# - all the end tags are at the end of the row with the information and not in the next row


# Example
#########
# srun --x11 --mem=8G --pty bash
# cd ~/work/hal.automate
# pgm=~/fragencode/tools/multi/Scripts/zotero2hal.awk
# time awk -v base=file -v fileRef1=irsd.authors.tsv -v fileRef2=template.art.xml -f $pgm zotero.output.xml > zotero.output.txt
# will produce 2 xml files that can then be sent with the curl command above
# file.1.xml and file.2.xml
# to see if they are correct just use the curl command above
# curl -X POST -d @file.1.xml -v -u ${HAL_USERNAME}:${HAL_PASSWORD} https://api.archives-ouvertes.fr/sword/hal/ -H "Packaging:http://purl.org/net/sword-types/AOfr" -H "Content-Type:text/xml"

# fileRef1=irsd.authors.tsv
# christophe	bureau
# jean		monlong
# 12 (2 fields)
# 1 (3 fields) *** between aguilar and martinez there is a space but all fornames and surnames are separated by a tab

# fileRef2=template.art.xml
# <?xml version="1.0" encoding="UTF-8"?>
# <TEI xmlns="http://www.tei-c.org/ns/1.0" xmlns:hal="http://hal.archives-ouvertes.fr/">
# 37 (1 fields)
# 11 (2 fields)
# 6 (3 fields)
# 3 (4 fields)
# 1 (5 fields)
# 1 (6 fields)
# 1 (11 fields)

# main input zotero.output.xml
# <?xml version="1.0" encoding="UTF-8"?>
# 
# 1 (0 fields)
# 493 (1 fields)
# 41 (2 fields)
# 4 (3 fields)
# 3 (4 fields)
# 2 (14 fields)

# main output file zotero.output.txt
# -             Predicting survival in patients with ‘non-high-risk’ acute variceal bleeding receiving β-blockers+ligation to prevent re-bleeding,               Balcar L,               Mandorfer M,               Hernández-Gea V,               Procopet B,               Meyer E,               Giráldez Á,               Amitrano L,               Villanueva C,               Thabut D,               Samaniego L,               Silva-Junior G,               Martinez J,               Genescà J,               Bureau C,               Trebicka J,               Herrera E,               Laleman W,               Palazón Azorín J,               Alonso J,               Gluud L,               Ferreira C,               Cañete N,               Rodríguez M,               Ferlitsch A,               Mundi J,               Grønbæk H,               Hernandez Guerra M,               Sassatelli R,               Dell'Era A,               Senzolo M,               Abraldes J,               Romero-Gómez M,               Zipprich A,               Casas M,               Masnou H,               Primignani M,               Krag A,               Nevens F,               Calleja J,               Jansen C,               Catalina M,               Albillos A,               Rudler M,               Tapias E,               Guardascione M,               Tantau M,               Schwarzer R,               Reiberger T,               Laursen S,               Lopez-Gomez M,               Cachero A,               Ferrarese A,               Ripoll C,               La Mura V,               Bosch J,               García-Pagán J,              Journal of Hepatology,               2024,             01688278,               80,               1,               73-81,             10.1016/j.jhep.2023.10.007
# -             Incidence and factors predictive of recurrent thrombosis in people with non-cirrhotic portal vein thrombosis,               Baiges A,               Procopet B,               Silva-Junior G,               Llop E,               Tellez L,               Darnell A,               Garcia-Criado Á,               Turon F,               Nicoara-Farcau O,               González-Alayón C,               Larrue H,               Magaz M,               Olivas P,               Perez-Campuzano V,               Calleja J,               Albillos A,               Reverter J,               Bureau C,               Bosch J,               Hernández-Gea V,               Garcia-Pagán J,              Journal of Hepatology,               2023,             01688278,               78,               1,               114-122,             10.1016/j.jhep.2022.08.023
# 1 (66 fields)
# 1 (139 fields)

# one of the TEI xml output file, eg file.1.xml
# <?xml version="1.0" encoding="UTF-8"?>
# <TEI xmlns="http://www.tei-c.org/ns/1.0" xmlns:hal="http://hal.archives-ouvertes.fr/">
# 210 (1 fields)
# 121 (2 fields)
# 51 (3 fields)
# 19 (4 fields)
# 2 (5 fields)
# 1 (6 fields)
# 1 (11 fields)
# 2 (16 fields)


# In the begin part we will:
############################
# - first read the list of irsd authors and store them in a table
# - second read the template xml file and remember each row (only when we will print the final files will
#   we look at the content of the templace rows to see whether we should substitute something with content
#   ({JOURNAL}, {DOI}, {AUTHORS}, {TITLE})
BEGIN{
    while (getline < fileRef1 >0)
    {
	l++;
	split($0,a,"\t");
	forename[l]=a[1];
	surname[l]=a[2];
    }
    # here l will be the total number of irsd authors we want to affiliate
    
    while (getline < fileRef2 >0)
    {
	m++;
	row[m]=$0;
    }
    # here m will be the number of rows in the template file
}


# The main script will read the main zotero xml file and will remember for each of them
# the minimal information that we need to put to the template file, each of them in a different
# table, with an index that corresponds to the rank of the article in the zotero file
# - title table with a single index, the rank of the article in the zotero file
# - journal table with a single index, the rank of the article in the zotero file
# - the doi table with a single index, the rank of the article in the zotero file
# - the authors table with a double index, one for the rank of the article in the zotero file and another one
#   for the rank of the author in the article in question
{
    # each time we see <biblStruct we are in a new article and therefore increment n
    if($0~/<biblStruct/)
    {
	n++;
    }
    # first we look for the title of the article and we know it is in the line that follows <title level="a">
    # so if we see this we assign okt to 1
    if(okt==0&&$0~/<title level/&&$0~/a/&&$0!~/j/)
    {
	okt=1;
    }
    # here we want the line just after <title level="a">
    if(okt==1&&$0!~/<title level/)
    {
	split($0,a,"<");
	title[n]=a[1];
	okt=0;
    }
    # second we look for the doi of the article and we know it is in the line that follows <idno type="DOI">
    # so if we see this we assign okd to 1
    if(okd==0&&$0~/<idno type/&&$0~/DOI/)
    {
	okd=1;
    }
    # here we want the line just after <idno type="DOI">
    if(okd==1&&$0!~/<idno/)
    {
	split($0,a,"<");
	doi[n]=a[1];
	okd=0;
    }
    # third we look for the journal information = title, idno, volume, issue, page and date
    # in variables called titjourn[n],idnojourn[n],voljourn[n],issjourn[n],pagejourn[n],datejourn[n]
    if(okj==0&&$0~/<title level/&&$0~/j/&&$0!~/a/)
    {
	okj=1;
    }
    if(okj==1&&$0!~/<title level/)
    {
	if($0~/title/)
	{
	    split($0,a,"<");
	    titjourn[n]=a[1];
	}
	else
	{
	    if(okj==1&&$0~/idno type/)
	    {
		lineidno=NR+1;
	    }
	    else
	    {
		if(okj==1&&NR==lineidno)
		{
		    split($0,a,"<");
		    idnojourn[n]=a[1];
		}
		else
		{
		    if(okj==1&&$0~/volume/)
		    {
			linevol=NR+1;
		    }
		    else
		    {
			if(okj==1&&NR==linevol)
			{
			    split($0,a,"<");
			    voljourn[n]=a[1];
			}
			else
			{
			    if(okj==1&&$0~/issue/)
			    {
				lineiss=NR+1;
			    }
			    else
			    {
				if(okj==1&&NR==lineiss)
				{
				split($0,a,"<");
				issjourn[n]=a[1];	
				}
				else
				{
				    if(okj==1&&$0~/page/)
				    {
					linepage=NR+1;
				    }
				    else
				    {
					if(okj==1&&NR==linepage)
					{
					    split($0,a,"<");
					    pagejourn[n]=a[1];
					}
					else
					{
					    if(okj==1&&$0~/<date>/)
					    {
						linedate=NR+1;
					    }
					    else
					    {
						if(okj==1&&NR==linedate)
						{
						    split($0,a,"<");
						    datejourn[n]=a[1]
						}
						else
						{
						    if($0!~/<imprint>/)
						    {
							okj=0;
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    # fourth we look for the forename and surname of each author of the current article (they go from 1 to nbauth[n])
    # in variables called forenameauth[n,k],surnameauth[n,k] and increment nbauth[n] at the same time
    # the analytic block is the one including the authors
    if($0~/<analytic>/)
    {
	oka=1;
    }
    if(oka==1&&$0~/<forename>/)
    {
	linefore=NR+1;
    }
    else
    {
	if(oka==1&&NR==linefore)
	{
	    nbauth[n]++;
	    split($0,a,"<");
	    forenameauth[n,nbauth[n]]=a[1];
	}
	else
	{
	    if(oka==1&&$0~/<surname>/)
	    {
		linesur=NR+1;
	    }
	    else
	    {
		if(oka==1&&NR==linesur)
		{
		    split($0,a,"<");
		    surnameauth[n,nbauth[n]]=a[1];
		}
		else
		{
		    # when we see </analytic> it means we are at the end of the author list
		    if(oka==1&&$0~/<\/analytic>/)
		    {
			oka=0;
		    }
		}
	    }
	}
    }
    
}
# at the end of reading the zotero file, n is the total number of articles
# this is the number of xml files that we will produce in the end part, see below

# In the end part we will write for each article of the zotero file a different xml file
# which is like the template but with the information from the paper in question at the right place
END{
    # for each article, numbered i (going from 1 to n)
    for(i=1; i<=n; i++)
    {
	# first we write the main output which is a list of all papers formatted as
	# - title, comma separated author list (author as surname fornamefirstletter), journal, year, idno, volume, issue, page, doi
	sauth="";
	for(k=1; k<=nbauth[i]; k++)
	{
	    # find first non space character in forename
	    n3=split(forenameauth[i,k],e,"");
	    q=1;
	    found3=0
	    while(found3==0&&q<=n3)
	    {
		if(e[q]!=" ")
		{
		    found3=1;
		}
		q++;
	    }
	    sauth=(sauth)(surnameauth[i,k]" "(e[q-1]))(", ");
	}
	print "- "title[i]", "sauth" "titjourn[i]", "datejourn[i]", "idnojourn[i]", "voljourn[i]", "issjourn[i]", "pagejourn[i]", "doi[i]; 
	# then we write the TEI hal compatible xml file for each paper, named after the article number in the zotero file (i here)
	f=base"."i".xml";
	# for each row j of the template from 1 to m, write it except when we need to replace something by a string (title, doi),
	# a block (journal) or several blocks (authors)
	for(j=1; j<=m; j++)
	{
	    r=row[j];
	    if(r~/TITLE/) # in this case we need to replace {TITLE} by the actual title (will be done twice for each article)
	    {
		split(r,a,"{TITLE}");
		print (a[1])(title[i])(a[2]) > f;
	    }
	    else
	    {
		if(r~/DOI/)  # in this case we need to replace {DOI} by the actual doi (will be done once for each article)
		{
		    split(r,a,"{DOI}");
		    print (a[1])(doi[i])(a[2]) > f;
		}
		else
		{
		    if(r~/JOURNAL/) # in this case we need to replace {JOURNAL} by the journal block so just write with the writejournal function
			            # will be done once for each article
		    {
			writejournal(titjourn[i],idnojourn[i],voljourn[i],issjourn[i],pagejourn[i],datejourn[i],f); 
		    }
		    else
		    {
			if(r~/AUTHORS/)  # in this case we need to replace {AUTHOR} by as many author blocks as there are authors so just write
			                 # with the writeauthor function as many times as there are authors
			                 # will be done once for each article
			{
			    # for each author from 1 to nbauth[i] for article i
			    for(k=1; k<=nbauth[i]; k++)
			    {
				writeauthor(forenameauth[i,k],surnameauth[i,k],f);
			    }
			}
			else	# in case there is nothing to replace in the template row, then we simply write the template row as such
			{
			    print r > f;
			}
		    }
		}
	    }
	}
    }
}

# This function takes as input 6 parameters about the journal of the article and 1 about the name of the file
# to write in, and writes the journal block in this file
function writejournal(tit,idno,vol,iss,pag,dat,fart){
    print "<idno type=\"issn\">"(idno)"</idno>" > fart;
    print "<title level=\"j\">"(tit)"</title>" > fart;
    print "<imprint>" > fart;
    print "  <biblScope unit=\"volume\">" > fart;
    print "  "(vol)"</biblScope>" > fart;
    print "  <biblScope unit=\"issue\">" > fart;
    print "  "(iss)"</biblScope>" > fart;
    print "  <biblScope unit=\"pp\">" > fart;
    print "  "(pag)"</biblScope>" > fart;
    print "  <date>" > fart;
    print "  "(dat)"</date>" > fart;
    print "</imprint>" > fart;
}

# This function takes as input 2 parameters about the author of an article, the forename and the surname
# and 1 about the name of the file to wirte in, and writes the author block in this file
# adding the IRSD affiliation in case the author belongs to it (provided in fileRef2)
function writeauthor(foren,surn,fart){
    print "<author role=\"aut\">" > fart;
    print "  <persName>" > fart;
    print "    <forename type=\"first\">"(foren)"</forename>" > fart;
    print "    <surname>"(surn)"</surname>" > fart;
    print "  </persName>" > fart;
    fo=tolower(foren);
    su=tolower(surn);
    found=0;
    p=1;
    while(found==0&&p<=l)
    {
	ifo=tolower(forename[p]);
	isu=tolower(surname[p]);
	n1=split(fo,c,ifo);
	n2=split(su,d,isu);
	if(n1>=2&&n2>=2)
	{
	    found=1;
	}
	p++;
    }
    # if(fo~/christophe/&&su~/bureau/)
    if(found==1)
    {
	print "  <affiliation ref=\"#struct-1002454\"/>" > fart;
    }
    print "</author>" > fart;
}
