#!/bin/sh

if [ -z "$1" ] ; then 
  echo "ERROR: Give the name of the CDS file as an argument";
  echo "The CDS file should be a gff file modified with the UCE ID as column one and the protein ID as column 10."
  echo "(columns are delimited by tabs)"
  echo "like:"
  echo 'uce-500	scaffold_10	JGI	CDS	526047	526419	.	+	0	exonNumber "1"; name "fgenesh2_pm.10_#_144"; proteinId "153259"	153259
'
  echo ""
  exit 1
fi

echo "Creating sqlite db with CDS info"
rm -f CDS.sqlite

cat << EOF >import.sql
DROP TABLE IF EXISTS CDS;
create table CDS (locus text, scaffold text, source text, type text, start int, end int, unk text, strand text, unk2 text, description text, proteinid text);
.separator "\t"
.import $1 CDS
EOF

sqlite3 CDS.sqlite <import.sql || echo "There was an error creating the DB of CDS"



echo "Exporting all the loci into the text file all_loci_CDS"

cat << EOF >export-all.sql
.mode tabs
.output all_loci_CDS
SELECT DISTINCT locus FROM CDS;
EOF

sqlite3 CDS.sqlite <export-all.sql || echo "There was an error exporting from the database"



echo "Finding loci where >1 proteinid was found"
while read locus; do
  count=$(sqlite3 CDS.sqlite "SELECT DISTINCT locus, proteinid from CDS WHERE locus = '$locus';" | wc -l)
  if [ $count -ne 1 ]; then
    echo "$locus: $count ProteinIDs found, removing this locus from the database"
    sqlite3 CDS.sqlite "DELETE FROM CDS WHERE locus = '$locus';"
  fi
done <all_loci_CDS

echo "Exporting UNIQUE loci into the text file loci_protein"

cat << EOF >export-unique.sql
.mode tabs
.output loci_protein
SELECT DISTINCT locus, proteinid FROM CDS;
EOF

sqlite3 CDS.sqlite <export-unique.sql || echo "There was an error exporting from the database"


echo "Cleaning up intermediate files (import.sql, export-all.sql, all_loci_CDS, export-unique.sql)"
rm import.sql export-all.sql all_loci_CDS export-unique.sql

echo "Output: loci_protein (Lists each UCE locus with one protein ID found, and its protein ID)"
