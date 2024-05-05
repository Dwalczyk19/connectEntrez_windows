#main
import sys
import ncbi.datasets.openapi    
from ncbi.datasets.openapi.api import gene_api
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi
from Bio import Entrez, SeqIO
from numpy import array as arr
import pandas as pd 
import os
import shutil
import random
import getpass
from zipfile import ZipFile
import jsonlines

class get_ids:       
    def download(self, zip_name, limits, df, m):
        if m == 1: 
            with DatasetsApiClient() as api_client:
                
                gene_api = DatasetsGeneApi(api_client)
                try:
                    print("Starting download of...{}".format(zip_name))
                    gene_dataset_download = gene_api.download_gene_package(
                        df[int(limits[0]):int(limits[1])],
                        include_annotation_type=["FASTA_CDS"],
                        _preload_content=False,
                    )
                    
                    with open(zip_name, "wb") as f:
                        f.write(gene_dataset_download.data)
                    print("Finished Downloading...")
                except DatasetsApiException as e:
                    sys.exit(f"Exception when calling GeneApi: {e}\n")
        elif m == 2:
            with DatasetsApiClient() as api_client:
                
                gene_api = DatasetsGeneApi(api_client)
                try:
                    print("Starting download of...{}".format(zip_name))
                    gene_dataset_download = gene_api.download_gene_package(
                        df,
                        include_annotation_type=["FASTA_CDS"],
                        _preload_content=False,
                    )
                    
                    with open(zip_name, "wb") as f:
                        f.write(gene_dataset_download.data)
                    print("Finished Downloading...")
                except DatasetsApiException as e:
                    sys.exit(f"Exception when calling GeneApi: {e}\n")
        elif m == 4: 
            with DatasetsApiClient() as api_client:
                
                
                gene_api = DatasetsGeneApi(api_client)
                gene_dataset_download = gene_api.download_gene_package(
                    df,
                    include_annotation_type=["FASTA_GENE", "FASTA_RNA"],
                    _preload_content=False,
                )
                
                with open(zip_name, "wb") as f:
                    f.write(gene_dataset_download.data)
                
        
            
    
    def get_gene(self,df): #df comes in as 1-d pandas df
        n = len(df)
        L = []
        if n > 1000: 
            if n % 1000 == 0: 
                steps =int( n / 1000)
                for i in range(1,steps+1):
                    ct = i * 1000
                    L.append((ct-1000,ct))
                
            else: 
                steps = int(n / 1000)
                rm = n % 1000
                
                if steps == 1: 
                    L.append((0,1000))
                    L.append((1000,1000+rm))
                else: 
                    for i in range(1,steps+1): 
                        ct = i * 1000 
                        L.append((ct-1000,ct))
                        if i == steps: 
                            L.append((ct, ct + rm))
                            
            
            for items in L: 
                num = f'{items[0], items[1]}'
                zipfile_name = "gene_cds" + num + ".zip"
                self.download(zipfile_name, items, df, 1)
            
                        
        else: 
            num = "(0, " + str(n) + ")"
            zipfile_name = "gene_cds" + num + ".zip"
            self.download(zipfile_name, None, df, 2)
        
    def checkAssembly(self,file, report_file):

        gene_file = {} 
    
        #file
        with jsonlines.open(report_file) as reader:
    
        
            for obj in reader.iter(type=dict):     
    
                try:
                    gene_file[obj["geneId"]] = [ obj["transcripts"][0]["cds"]["accessionVersion"]] #length one 
                    
                except KeyError as e: 
                    e = str(e)
                    
                    if e == "'cds'":
                        if obj["geneId"] == "55199": #change to iterate over all types of genes (this is psuedogene & id specific)
                            gene_file[obj["geneId"]] = [obj["transcripts"][0]["type"], "FASTA_GENE"]
                            
                        else: 
                            gene_file[obj["geneId"]] = [obj["type"], "CDS"]
                        
                    
                    elif e == "'transcripts'":
                        gene_file[obj["geneId"]] = [obj["type"], "FASTA_GENE"] #length two with FASTA_GENE at index 1
            
            return gene_file
        
    def findAssembly(self, L, home):  #convert files in L to something that can be altered and changed into a sequence. Maybe make an initial step to check the validity of each gene via some datasets function
        #for .jsonl
        #home = dir_path = r'users\david\downloads\entrezPackage
        
        complete = pd.DataFrame()
        from operator import itemgetter
        final = pd.DataFrame()
        for i in range( len( L )): 
            print(L[i], "============================================")
            with ZipFile(L[i], "r") as zip: 
                direct = zip.namelist()
                report = direct[2]
                cds_file = direct[1]
                zip.extract(report, path = "assembly" + str(i))
                zip.extract(cds_file, path = "assembly" + str(i))
                
            new_file = home + "\\assembly" + str(i)+ "\\ncbi_dataset\\data\\data_report.jsonl"
            cds_file = home + "\\assembly" + str(i) + "\\ncbi_dataset\\data\\" + "cds.fna"
            gene_file = self.checkAssembly(L[i], new_file)
            genesort = list(gene_file.keys())
            genesort.sort()
            gene_dict = {int(i): gene_file[i] for i in genesort}
            cds = [y[0] for x,y in gene_dict.items() if len(y) == 1]
            assembly = [x for x,y in gene_dict.items() if (len(y) > 1 and y[1] == "CDS")]
            gene_fna = [x for x,y in gene_dict.items() if (len(y) > 1 and y[1] == "FASTA_GENE")]
            cds_df = pd.DataFrame()
            test_df = pd.DataFrame()
            gene_df = pd.DataFrame()
            with open(cds_file, "r") as handle: #main part
                final = []
                check = 0
                count = 0
                find_length = []
                for record in SeqIO.parse(handle, "fasta"):
                    
                    transcript = str(record.id).split(":")[0] #NM/XM, NR/XR
                    geneID = str(record.description).split("] [")[1][7:] #Gene ID
                    geneName = str(record.description).split()[1].strip() #Gene name 
                    seq = str(record.seq) #sequence
                    length = len(seq) #length
                    
                    if transcript in cds: 
                        cds_df = pd.concat([cds_df, pd.DataFrame({"Transcript": [transcript], "ID": geneID, "Name":geneName, "Sequence": seq, "Length": length})], ignore_index=True)
                    
                    
                    elif int(geneID) in assembly:
                        check = int(geneID)
                        if len(find_length) == 0:
                            find_length.append(check)
                            find_length.append( (transcript, length, geneName, seq))
                        elif check in find_length:
                            find_length.append((transcript, length, geneName, seq))
                        elif check not in find_length:
                            name = find_length[0]
                            find_length.pop(0)
                            top = max(find_length, key = itemgetter(1))
                            test_df = pd.concat([test_df, pd.DataFrame({"Transcript": [top[0]], "ID":name, "Name": top[2], "Sequence": top[3], "Length":top[1]})], ignore_index=True)
                            find_length = []
                    
                    #elif int(geneID) in gene_fna:
                    #    print(int(geneID))
                        
                merge = pd.concat([test_df, cds_df])
                complete = pd.concat([complete, merge])
        
        complete[["ID", "Name", "Sequence"]].to_csv("Human-Entrez-IDs.csv", index = False)
        return complete



if __name__ == "__main__": 
    key = "9cbb475748fcce2a676126874b4fc0616f08"
    dir_path = "C:\\Users\\" + str(getpass.getuser()) + "\\Downloads\\entrezPackage" #main directory path where gene_cds files can be accessed
    os.chdir(dir_path)
    
    path2file = input("Enter file path: ")
    #path2file = "C:\\Users\\david\\entrezPackage\\Human_genes_with_entrez_IDs060523.xlsx"
    #read csv vs.xlsx
    if path2file.split(".")[1] == "csv": 
        df = pd.read_csv(path2file) 
    elif path2file.split(".")[1] == "xlsx":
        df = pd.read_excel(path2file)

    #if # of cols is > 1 
    if df.shape[1] > 1: 
        df = df.iloc[:,0]
    
    #apply class
    get = get_ids()
    get.get_gene(pd.unique(df).tolist())
    complete = os.listdir(dir_path)
    idx = [complete[complete.index(items)]  for items in complete if "gene_cds" in items]
    complete_df = get.findAssembly(idx, dir_path)
    
    #remove all useless files
    #doesn't work w/out admin privileges
    '''
    for items in os.listdir(dir_path):
        if "gene_cds" in items or "assembly" in items: 
            os.remove(dir_path + "\\"+ items)
    '''
    