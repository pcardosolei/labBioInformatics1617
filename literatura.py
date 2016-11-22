from Bio import Entrez
from Bio import Medline

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


def recordIdList(user,termo,num):
  idList = obterListaTermo(user,termo,num)
  if idList == -1:
    print("Não existem artigos")
  else:
    print(idList)

def recordMedline(user,termo,num):
  idList = obterListaTermo(user,termo,num)
  if idList == -1:
    print("Não existem artigos")
  else:
    handle = Entrez.efetch(db="pubmed",id=idList,rettype="medline",retmode="text", usehistory="y")
    records = Medline.parse(handle)
    # transformar isto numa lista
    records = list(records)
    print(records)

def test():
  recordMedline("history.a61043@alunos.uminho.pt","Legionella gene",10)

test()

"""
  VER a Literatura associada com o mesmo
"""
