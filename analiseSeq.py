"""
  Parte BioPython
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import ExPASy
from Bio import SwissProt


"""
  Imports genericos
"""
import os
import re
import time
import xml.etree.ElementTree as ET


class Feature:
  def __init__(self,ref):
    self.ref = ref

  #adicionar a funcao na classe feature
  def addFunction(self,function):
    self.function = function

  def addNote(self,note):
    self.note = note

  #adicionar a locusTag na classe feature
  def addDbRef(self,db_xref):
    self.db_xref = db_xref

  #adicionar a sequencia na classe feature
  def addSeq(self,seq):
    self.seq = seq

  def addProteinID(self,ident):
    self.id = ident

  def addProduct(self,prod):
    self.product = prod

  def addRefSeq(self,ref):
    self.ref = ref
""" ###########

  Aspetos gerais na analise da sequencia e features presentes na NCBI
""" ###########


"""
  Funcao trabalha com a sequencia.
  Le as features e guarda as informacoes que sao consideradas importantes - db_ref , function e locus_tag
"""
def locateFeatures():
  records = SeqIO.read(open("example.gb"),"genbank")
  dic = {}
  flag = False
  for feature in records.features:
    #vale a pena considerar a source e o STS presente na nossa proteina
    if "gene" in feature.type and flag == False:
      qualifiers = feature.qualifiers
      flag = True
      if "locus_tag" in qualifiers:
        aux = Feature(qualifiers['locus_tag'][0])
        if "db_xref" in qualifiers.keys():
          aux.addDbRef(qualifiers['db_xref'][0])
        tag = qualifiers['locus_tag'][0]
        dic[tag] = aux
    elif flag == True:
      qualifiers = feature.qualifiers
      flag = False
      if "locus_tag" in qualifiers:
        aux = qualifiers['locus_tag'][0]
        if "function" in qualifiers.keys():
          dic[aux].addFunction(qualifiers['function'])
        if 'translation' in qualifiers.keys():
          dic[aux].addSeq(qualifiers['translation'][0])
        if 'note' in qualifiers.keys():
          dic[aux].addNote(qualifiers['note'][0])
        if 'protein_id' in qualifiers.keys():
          dic[aux].addProteinID(qualifiers['protein_id'][0])
        if 'product' in qualifiers:
          dic[aux].addProduct(qualifiers['product'][0])
    else:
      pass
  return dic

"""
  read XML from blast
"""
def readXMLBlast(file):
  root = xml.etree.ElementTree.parse('blasts/GeneID:19834107.xml').getroot()
  print(root.findall("./BlastOutput_iterations"))


"""
  Funcao para obter valores adicionais de cada um dos genes -- ja nao me lembro onde xD
"""
def getInfoFeature(feature,tag):
  user  = "history.a61043@alunos.uminho.pt"
  if hasattr(feature,"id"):
    pass
    #writeProteinTag(user,feature.id,tag)
  else:
    #ESTA PARTE DO CODIGO NAO ESTA A SAKAR PARTES CORRECTAS - SAKA GENES DIFERENTES
    #COLOCAR UM PRINT E SABER QUAIS SAO E USAR A FUNCAO INDIVIDUAL PARA SAKAR OS GENES QUE FALTAM WriteInfoTag
    genID = re.split(r':',feature.db_xref)
    writeProteinTag(user,genID,tag)



"""
  UNIPROT Functions
"""
def getUniProt():
  """
    obter informacao uniprot de uma unica prot
  """
  pass


def getUniProtLegionella(dic):
  """
    Procurar tdas as proteinas da legionella na uniprot
  """
  for tag in dic:
    getUniProt(dic[tag].db_xref)



"""
  funcao para realizar os blasts as sequencias dos locustag
"""
def blastAllSeq(folder):
  soma = 0
  fastas = []
  fastas += [each for each in os.listdir(folder) if each.endswith('.fasta')]

  for fasta in fastas:

    fastaToBlastFile(fasta)
    soma += 1
    print("a Descansar")
    time.sleep(180) #pausar a funcao por 5 minutos para nao estourar
    print("Proteina {}".format(soma))
  print("Foram tratadas {} sequencias".format(soma))




"""
  Passagem de informacao das features
"""
def stringInfo(feature):
  text = ""
  if hasattr(feature,"db_xref"):
    text += "db_xref: "+ feature.db_xref + "\n"
  if hasattr(feature,"function"):
    text += "Function: "+ ''.join(feature.function) + "\n"
  if hasattr(feature,"note"):
    text += "Note: " + feature.note + "\n"
  if hasattr(feature,"id"):
    text += "Protein_id: " + feature.id + "\n"
  if hasattr(feature,"product"):
    text += "Product: " + feature.product + "\n"
  text += " \n"
  return text

"""
  Comparar a ProteinTable da legionella com os resultados de todas as proteinas
  PODE
"""
def listProteinTable(dicto):
  with open("ProteinTable416_0.txt") as fp:
    flag = 0
    ofile = open("teste.txt","w")
    for line in fp:
      if flag > 0:
        text = line.split()
        aux = text[11:]
        genID = "GeneID:"+text[5]
        if genID in dicto:
          ofile.write(genID+"\n")
          ofile.write(' '.join(aux)+"\n")
          if hasattr(dicto[genID],'product'):
            func = ''.join(dicto[genID].product)
            ofile.write(func+"\n\n")
          else:
            ofile.write("\n")
      flag += 1



def test(user):
  entrezQuery = "refseq[filter] AND txid%s"%(ncbiTaxId)
  searchResultHandle = Entrez.esearch(db="protein", term=entrezQuery)
  searchResult = Entrez.read(searchResultHandle)
  searchResultHandle.close()

######################################################################################################################
"""
  ITERATORS
"""

def writeGeneXML(user,dic):
  """
    Iterador para ir buscar o xml correspondente de cada gene da legionella
  """
  for feature in dic:
    genID = re.split(r':',dic[feature].db_xref)
    writeSequenceXML(user,genID[1],feature)

def getProteinFeatures(dic):
  """
  Iterador de cada uma das proteinas na legionella
  """
  for feature in dic:
    getInfoFeature(dic[feature],feature)

######################################################################################################################
"""
  LEITURAS
"""

def readNCBI(user,startSeq,stopSeq):
  """
  Ler a sequencia do NCBI da legionella

  parametros - utilizador , posicao inicial da sequencia , posicao final da sequencia
  """
  Entrez.email = user
  handle = Entrez.efetch(db="nucleotide",id="NC_002942.5",rettype="gbwithparts",retmode="text",seq_start=startSeq,seq_stop=stopSeq)
  return handle

def readNCBIGene(user,idG):
  """
    Le informacao gene
  """
  Entrez.email = user
  handle = Entrez.efetch(db="gene",id=idG,rettype="gb",retmode="xml")
  return handle


def readProteinTag(user, genID):
  """
    funcao para ir buscar mais informacoes de cada uma das locus tags
  """
  Entrez.email = user
  handle = Entrez.efetch(db="Genome",id=genID,rettype="gbwithparts",retmode="text")
  return handle

######################################################################################################################


"""
  Funcoes de mexer com o xml
"""

def readXMLGene(folder,dic):
  xmls = []
  xmls += [each for each in os.listdir(folder) if each.endswith('.xml')]
  for xml in xmls:
    ofile = open(folder+"/"+xml,"r")
    tree = ET.parse(ofile)
    root = tree.getroot()
    place = root.findall("./Entrezgene/Entrezgene_comments/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_str")
    locustag = re.split(r'\.',xml)[0]
    for pl in place:
      dic[locustag].addRefSeq(pl.text)







######################################################################################################################


"""
  ESCRITAS
"""
def writeSequence(user,startSeq,stopSeq):
  """
   escrever a sequencia total da legionella para um ficheiro
  """
  handle = readNCBI(user,startSeq,stopSeq)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(handle,"teste.gb","xml")

def writeSequenceXML(user,idG,tag):
  """
   escrever a sequencia total da legionella para um ficheiro
  """
  handle = readNCBIGene(user,idG)
  ofile = open("geneXML/"+tag+".xml","w")
  ofile.write(handle.read())
  ofile.close()


def writeProteinTag(user,genID,tag):
  """
  escrever a informacao da locus tag para cada um dos ficheiros
  """
  handle = readProteinTag(user,genID)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(record,"locTags/"+tag+".gb","genbank")

def proteinSequencesToFasta(dic):
  """
  Funcao que pega nas proteinas que tem uma sequencia associada e cria um ficheiro FASTA
  """
  for tag in dic:
    if hasattr(dic[tag], 'seq'):
      ofile = open("fastas/"+tag+".fasta", "w")
      ofile.write(">" + tag + "\n" +dic[tag].seq + "\n")
      #do not forget to close it
      ofile.close()
    else:
      pass

def fastaToBlastFile(fasta):
  """
    Pega no blast da sequencia e vai buscar o blast associado e guarda o
  """
  ofile = open("fastas/"+fasta,"r")
  print("Reading .... {}".format(ofile))
  seq = ""
  for line in ofile:
    seq += line
  result_handle = NCBIWWW.qblast("blastp", "nr",seq)
  name = os.path.basename(fasta)
  blastFile = open("blasts/"+os.path.splitext(name)[0]+".xml","w")
  blastFile.write(result_handle.read())
  blastFile.close()
  result_handle.close()
  ofile.close()



def saveGenes(dictionary):
  """
    Guardar informacao genes
  """
  ofile = open("genes.txt","w")
  for gene in dictionary:
    txt = stringInfo(dictionary[gene])
    ofile.write(gene+"\n")
    ofile.write(txt)
  ofile.close()


######################################################################################################################


"""
  MISC
"""


def createFolders(directories):
  """
    Create folders
  """
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)



######################################################################################################################


"""
  MAIN
"""
def main():
  #readXMLBlast("1")

  #folders = ["fastas","blasts","locTags"]
  #createFolders(folders)
  #writeSequence("history.a61043@alunos.uminho.pt","2873801","3124550")
  dic = locateFeatures()
  #writeGeneXML("history.a61043@alunos.uminho.pt",dic)
  #getProteinFeatures(dic)
  #saveGenes(dic)
  #fastaToBlastFile("fastas/GeneID:19834107.fasta")
  #blastAllSeq("fastas")
  #getSwissProtInfo("19834166")
  #readBlastProtein()
  #getUniProt()
  #test("history.a61043@alunos.uminho.pt")

  readXMLGene("geneXML",dic)

"""
  TODO - read protein locus genbanks
       - read BLAST FILES
"""
main()
