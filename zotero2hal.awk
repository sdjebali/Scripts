# zotero2hal.awk
# the idea is to automate the submission of papers from our team (specially clinical papers
# with a lot of authors) into hal for the hceres
# This script starts from two inputs:
# - fileRef is a template xml file to submit a paper in hal, with 5 places that should be replaced with content
#   * {TITLE} at two places with a string each time (the same one)
#   * {AUTHORS} at one place with several blocks of authors
#   * {JOURNAL} at one place with a block
#   * {DOI} at one place with a string
# - main input file that is the zotero output in xml format with several papers that we want to submit to hal
# This script will produce as many hal compatible xml files as we have papers in the main input file
# The idea is then to submit all of them one by one using a command like this one
# curl -X POST -d @art1.works.xml -v -u ${HAL_USERNAME}:${HAL_PASSWORD} https://api.archives-ouvertes.fr/sword/hal/ -H "Packaging:http://purl.org/net/sword-types/AOfr" -H "Content-Type:text/xml"
# where art1.works.xml is one of the xml file produced by this script

# Several assumptions are done here:
####################################
# - the title of the article is in the row that follows the row <title level="a">
# - the doi of the article is in the row that follows the row <idno type="DOI">


# Example
#########
# cd ~/work/hal.automate
# pgm=~/fragencode/tools/multi/Scripts/zotero2hal.awk
# time awk -v fileRef=template.art.xml -f $pgm zotero.output.xml
# will produce 2 xml files that can then be sent with the curl command above


# in the begin part we will read the template xml file and remember each row
# as well as the rows where we should put information from the zotero file
BEGIN{
    while (getline < fileRef >0)
    {
	m++;
	row[m]=$0;
   }
}
# at the end here, m is the number of rows in the template file

# the main script will read the main zotero xml file and will remember for each of them
# the minimal information that we need to put to the template file, each of them in a different
# table, with an index that corresponds to the rank of the article in the zotero file
# - title table with a single index, the rank of the article in the zotero file
# - journal table with a single index, the rank of the article in the zotero file
# - the doi table with a single index, the rank of the article in the zotero file
# - the authors table with a double index, one for the rank of the article in the zotero file and another one
#   for the rank of the author in the article in question
{
    # each time we see <biblStruct we are in a new article and therefore increment n
    if($0~/\<biblStruct/)
    {
	n++;
    }
    # first we look for the title of the article and we know it is in the line that follows <title level="a">
    # so if we see this we assign okt to 1
    if($0~/\<title level/&&$0~/a/&&$0!~/j/)
    {
	okt=1;
    }
    # here we want the line just after <title level="a">
    if(okt==1&&$0!~/title/)
    {
	split($0,a,"<");
	title[n]=a[1];
	okt=0;
    }
    # second we look for the doi of the article and we know it is in the line that follows <idno type="DOI">
    # so if we see this we assign okd to 1
    if($0~/\<idno type/&&$0~/DOI/)
    {
	okd=1;
    }
    # here we want the line just after <idno type="DOI">
    if(okd==1&&$0!~/idno/)
    {
	split($0,a,"<");
	doi[n]=a[1];
	okd=0;
    }
    # third we look for the journal information = title, idno, volume, issue, page and date
    # in variables called titjourn[i],idnojourn[i],voljourn[i],issjourn[i],pagejourn[i],datejourn[i]

    # fourth we look for the forname and surname of each author
    # in variables called forenauth[k],surnameauth[k] for k from 1 to nbauth[i]
    
}
# at the end of reading the zotero file, n is the total number of articles
# this is the number of xml files that we will produce in the end part, see below

# In the end part we will write for each article of the zotero file a different xml file
# which is like the template but with the information from the paper in question at the right place
END{
    # for each article i from 1 to n
    for(i=1; i<=n; i++)
    {
	# a different xml file named after the article rank in the zotero file
	f="file."i".xml";
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
				writeauthor(forenauth[k],surnameauth[k],f);
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
    print "author role=\"aut\">" > fart;
    print "  <persName>" > fart;
    print "    <forename type=\"first\">"(foren)"</forename>" > fart;
    print "    <surname>"(surn)"</surname>" > fart;
    print "  </persName>" > fart;
    print "</author>" > fart;
}
