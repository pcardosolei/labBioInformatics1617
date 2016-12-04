from Bio import Entrez
from Bio import Medline

import os
import json

"""
  Esta funcao serve para encontrar os IDs do artigos que estão relacionados com a frase colocada
"""
def obterListaTermo(user,termo,num):
  Entrez.email = user
  handle = Entrez.egquery(term=termo, usehistory="y")
  record = Entrez.read(handle)

  """
    Relaxamento na busca de artigos
  """
  for row in record["eGQueryResult"]:
    if row["DbName"] == "pubmed":
      value = row["Count"]
      if value == 0:
        return -1
      elif int(value) > num:
        value = str(num)

  handle = Entrez.esearch(db="pubmed", term=termo,retmax=value, usehistory="y")
  record = Entrez.read(handle)
  idList = record["IdList"] #todos os ids com o termo legionella
  return idList

"""
  Funcao para obter informacao do artigo a partir do ids da funcao em cima
"""
def recordMedline(user,termo,num):
  idList = obterListaTermo(user,termo,num)
  if idList == -1:
    print("Não existem artigos")
    return -1
  else:
    handle = Entrez.efetch(db="pubmed",id=idList,rettype="medline",retmode="text", usehistory="y")
    records = Medline.parse(handle)
    # transformar isto numa lista
    records = list(records)
    print(records)
    return records

"""
  Gravar Ficheiro Literatura
"""
def saveLiteratura(termo,records):
  ofile = open("literatura/"+termo+".txt","w")
  ofile.write(termo+"\n\n")
  for record in records:
    ofile.write(json.dumps(record)+"\n\n")
  ofile.close()

"""
  Criar pasta Literatura
"""
def createFolder(directory):
    if not os.path.exists(directory):
      os.makedirs(directory)

"""
  Zona de teste
"""
def main():
  termo = "Legionella gene" #frase a usar nas pesquisas
  createFolder("literatura")
  records = recordMedline("history.a61043@alunos.uminho.pt",termo,10)
  saveLiteratura(termo,records)

########
main()


