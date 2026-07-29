from Bio import Entrez
import pandas as pd
import time

Entrez.email = "rackbersingh@outlook.com"  # Replace with your email

accessions = [
    "GCF_001274345.1", "GCF_002402265.1", "GCF_004787495.1", "GCF_000837865.1",
    "GCF_000863805.1", "GCA_031117335.1", "GCA_031121055.1", "GCF_000859985.2",
    "GCF_000838265.1", "GCF_000186165.1", "GCF_000191725.1", "GCF_000005845.2",
    "GCF_000007325.1", "GCF_000007785.1", "GCF_000008345.1", "GCF_000008525.1",
    "GCF_000009585.1", "GCF_000009925.1", "GCF_000010505.1", "GCF_000015785.2",
    "GCF_000021265.1", "GCF_000024105.1", "GCF_000024565.1", "GCF_000024945.1",
    "GCF_000027165.1", "GCF_000027305.1", "GCF_000092265.1", "GCF_000146345.1",
    "GCF_000153885.1", "GCF_000153925.1", "GCF_000154805.1", "GCF_000156595.1",
    "GCF_000156615.2", "GCF_000156675.1", "GCF_000157995.1", "GCF_000158235.1",
    "GCF_000159115.1", "GCF_000159695.1", "GCF_000159995.1", "GCF_000160035.2",
    "GCF_000160055.1", "GCF_000160075.2", "GCF_000160435.1", "GCF_000162015.1",
    "GCF_000163475.1", "GCF_000164675.2", "GCF_000164695.2", "GCF_000173355.1",
    "GCF_000173915.1", "GCF_000174015.1", "GCF_000174175.1", "GCF_000175635.1",
    "GCF_000177075.1", "GCF_000177375.1", "GCF_000183565.2", "GCF_000191725.2",
    "GCF_000195995.1", "GCF_000196035.1", "GCF_000196535.1", "GCF_000208405.1",
    "GCF_000212375.1", "GCF_000239695.1", "GCF_000253395.1", "GCF_000262545.1",
    "GCF_000296505.1", "GCF_000314675.2", "GCF_000315465.1", "GCF_000372405.1",
    "GCF_000374605.1", "GCF_000375645.1", "GCF_000376645.1", "GCF_000393015.1",
    "GCF_000411175.1", "GCF_000411375.1", "GCF_000412995.1", "GCF_000420065.1",
    "GCF_000420445.1", "GCF_000423485.1", "GCF_000430525.1", "GCF_000455445.1",
    "GCF_000463505.1", "GCF_000478985.1", "GCF_000510445.1", "GCF_000516535.1",
    "GCF_000565015.1", "GCF_000613345.1", "GCF_000724605.1", "GCF_000758765.1",
    "GCF_000818035.1", "GCF_001190755.1", "GCF_001298465.1", "GCF_001517935.1",
    "GCF_001552035.1", "GCF_001553115.1", "GCF_002214645.1", "GCF_002998925.1",
    "GCF_003019295.1", "GCF_006094375.1", "GCF_014202695.1", "GCF_014647755.1",
    "GCF_016127955.1", "GCF_019048645.1", "GCF_019931005.1", "GCF_020097375.1",
    "GCF_023614525.1", "GCF_030408675.1", "GCF_040556925.1", "GCF_900446675.1",
    "GCF_900453895.1", "GCF_900637555.1", "GCF_902374465.1", "GCF_943169685.1",
    "GCF_943169825.2", "GCF_943169845.1", "GCF_003033055.1", "GCF_000861845.1",
    "GCA_048568685.1", "GCF_000864765.1", "GCF_900090045.1", "GCF_000861085.1",
    "GCF_000864745.1", "GCF_000195955.2", "GCF_000410535.2"
]


ranks_of_interest = [
    "superkingdom", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "subspecies"
]

data = []

for acc in accessions:
    try:
        # Search for the assembly
        handle = Entrez.esearch(db="assembly", term=acc)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            assembly_id = record["IdList"][0]

            # Get assembly summary
            summary_handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(summary_handle)
            summary_handle.close()
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            taxid = docsum['Taxid']

            # Get taxonomy
            tax_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            tax_records = Entrez.read(tax_handle)
            tax_handle.close()
            tax_info = tax_records[0]

            # Build taxonomic lineage dictionary
            tax_lineage = {rank: "" for rank in ranks_of_interest}
            for item in tax_info.get("LineageEx", []):
                rank = item["Rank"]
                name = item["ScientificName"]
                if rank in tax_lineage:
                    tax_lineage[rank] = name

            tax_dict = {
                "Accession": acc,
                "AssemblyName": docsum.get("AssemblyName", ""),
                "AssemblyStatus": docsum.get("AssemblyStatus", ""),
                "AssemblySpan": docsum.get("AssemblySpan", ""),
                "ChromosomeCount": docsum.get("ChromosomeCount", ""),
                "ScaffoldCount": docsum.get("ScaffoldCount", ""),
                "ContigCount": docsum.get("ContigCount", ""),
                "RefSeqCategory": docsum.get("RefSeqCategory", ""),
                "WGSProject": docsum.get("WGSProject", ""),
                "Submitter": docsum.get("Submitter", ""),
                "SubmissionDate": docsum.get("SubmissionDate", ""),
                "Organism": docsum.get("Organism", ""),
                "SpeciesName": docsum.get("SpeciesName", ""),
                **tax_lineage,
                "ScientificName": tax_info.get("ScientificName", ""),
                "TaxID": tax_info.get("TaxId", ""),
                "Rank": tax_info.get("Rank", ""),
                "Division": tax_info.get("Division", ""),
                "GeneticCode": tax_info.get("GeneticCode", {}).get("GCName", ""),
                "MitoGeneticCode": tax_info.get("MitoGeneticCode", {}).get("MGCName", ""),
                "CommonNames": "; ".join(
                    tax_info.get("OtherNames", {}).get("GenbankCommonName", "").split(",")
                )
            }

            data.append(tax_dict)
            time.sleep(0.4)  # Be kind to NCBI's servers

    except Exception as e:
        print(f"Failed to fetch data for {acc}: {e}")
        continue

# Save to Excel
df = pd.DataFrame(data)
df.to_excel("detailed_virus_taxonomy.xlsx", index=False)
print("Saved detailed genome taxonomy info to 'detailed_virus_taxonomy.xlsx'")
