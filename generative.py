import random
import pickle
import pandas as pd
import json

import warnings
warnings.filterwarnings('ignore')

try:
    import docking
    df = pd.read_csv('drugbank_drug_target_label_mapping_amino_acid_pair_manage2_vers.csv')
except:
    import AIGDrugDeNovo.docking as docking
    # if in google colab environment
    df = pd.read_csv('/content/AIGDrugDeNovo/drugbank_drug_target_label_mapping_amino_acid_pair_manage2_vers.csv')
df['TargetName'] = df['TargetName'].dropna()
df = df[df['DrugSeq'] != '-----'][df['Drug_Class_1'] != 185]
df = df.reset_index(drop=True)
try:
    with open("drug_class_mapping.json", "r") as file:
        dc = json.load(file)
except:
    # in google colab envinronment
    filename = "/content/AIGDrugDeNovo/drug_class_mapping.json"
    with open(filename, "r") as file:
        dc = json.load(file)

dc = {str(k).lower():v for k,v in dc.items()}
dcr = {v:k for k,v in dc.items()}
#print(dcr)

def count_alphabets(string):
  """
  Menghitung jumlah alphabet "ACDEFGHIKLMNPQRSTVWY" dari sebuah input string.

  Args:
    string: String yang akan dihitung jumlah alphabetnya.

  Returns:
    Dictionary jumlah huruf dari string.
  """

  alphabets = "ACDEFGHIKLMNPQRSTVWY"
  counts = {}
  for char in string:
    if char in alphabets:
      counts[char] = counts.get(char, 0) + 1

  alphabets_ct = list(counts.keys())
  for char in alphabets:
    if char not in alphabets_ct:
      counts[char] = 0
  return counts

def exctraction_feature(string):
    feature = count_alphabets(string)
    feature = [feature[aa] for aa in "ACDEFGHIKLMNPQRSTVWY"]
    return feature

def predict(string):
    try:
      filename = "drug_nb_model.pkl"

      # Buka file pickle
      with open(filename, "rb") as file:
          # Load model dari file
          model = pickle.load(file)

      res = model.predict([exctraction_feature(string)])
    except:
      filename = "/content/AIGDrugDeNovo/drug_nb_model.pkl"

      # Buka file pickle
      with open(filename, "rb") as file:
          # Load model dari file
          model = pickle.load(file)

      res = model.predict([exctraction_feature(string)])

    return res

def shuffle_alphabet(string):
  """
  Mengacak posisi alphabet dari sebuah string.

  Args:
    string: String yang ingin diacak.

  Returns:
    String dengan posisi alphabet yang telah diacak.
  """

  alphabet = list(string)
  random.shuffle(alphabet)
  return ''.join(alphabet)

def is_equal_unordered(list1, list2, tolerance_num=1):
  """
  Fungsi untuk membanding 2 list yang sama tapi tidak berurutan

  Args:
    list1: List pertama
    list2: List kedua

  Returns:
    True jika list1 dan list2 sama, False jika tidak
  """

  # Ubah list menjadi set untuk menghilangkan urutan
  set1 = set(list1)
  set2 = set(list2)

  # Hitung jumlah value yang sama
  intersection_size = len(set1 & set2)
  # print(intersection_size)

  # Kembalikan True jika jumlah value yang sama minimal 4
  return intersection_size >= tolerance_num


class Generative:
  def __init__(self, drug_class=[], len_result=1):
     self.drug_class = drug_class
     self.len_result = len_result

  def run(self):
    print(self.drug_class)
    result = []
    result_target = {}
    dock_res = {}
    while(True):
      index_rd = random.randint(0, len(df)-1)
      seq = str(df['DrugSeq'][index_rd])
      print(seq)
      seq = shuffle_alphabet(seq)
      pred = predict(seq)
      pred = pred[0]
      # print(is_equal_unordered(self.drug_class , pred, 2))

      if is_equal_unordered(self.drug_class , pred, 1):
          pred_res = [dcr[i] for i in pred]
          result.append(seq)
          target_seq = str(df['TargetSequenceAminoAcid'][index_rd])
          result_target[str(df['TargetName'][index_rd])] = target_seq
          dock_res[str(df['TargetName'][index_rd])] = docking.docking(seq, target_seq).run()

      if len(result) == self.len_result:
          break

    return result, result_target, dock_res
