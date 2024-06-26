try:
    import generative
except:
    import AIGDrugDeNovo.generative
import sys
import json
import datetime
import random
import string

def generate_random_code():
    """
    Generates a random code with the format "aidrugdenovo_{date}_{time}_{random_code}".
    
    :return: A string representing the random code.
    """
    # Get current date and time
    now = datetime.datetime.now()
    date_str = now.strftime("%Y%m%d")
    time_str = now.strftime("%H%M%S")
    
    # Generate a random alphanumeric code
    random_code = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    
    # Format the code
    random_code_str = f"aidrugdenovo_{date_str}_{time_str}_{random_code}"
    
    return random_code_str
  
def export_dict_to_json(data, filename):
    """
    Exports a dictionary to a JSON file.

    :param data: The dictionary to export.
    :param filename: The name of the file to export to.
    """
    try:
        with open(filename, 'w') as json_file:
            json.dump(data, json_file, indent=4)
        print(f"Dictionary successfully exported to {filename}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Definisikan nama file json
try:
    filename = "drug_class_mapping.json"

    # Buka file json
    with open(filename, "r") as file:
        drug_class_dict = json.load(file)
except:
    # in google colab
    filename = "/content/AIGDrugDeNovo/drug_class_mapping.json"

    # Buka file json
    with open(filename, "r") as file:
        drug_class_dict = json.load(file)

drug_class_dict = {str(x).lower():y for x,y in drug_class_dict.items()}

# drug_class_dict = {str(x).lower():y for x,y in drug_class_dict.items()}
# # print(drug_class_dict)

# arguments = sys.argv[1:]
# dc1 = "none"
# dc2 = "none"
# dc3 = "none"
# dc4 = "none"
# dc5 = "none"
# dc6 = "none"
# dc7 = "none"

# for arg in arguments:
#   key, value = arg.split("=")
#   if key == "dc1":
#     dc1 = value
#   elif key == "dc2":
#     dc2 = value
#   elif key == "dc3":
#     dc3 = value
#   elif key == "dc4":
#     dc4 = value
#   elif key == "dc5":
#     dc5 = value
#   elif key == "dc6":
#     dc6 = value
#   elif key == "dc7":
#     dc7 = value
#   else:
#     print("Perintah tidak dikenali")

# drug_class = [dc1, dc2, dc3, dc4, dc5, dc6, dc7]
# drug_class = [drug_class_dict[dc] for dc in drug_class]
# print(drug_class)
# generate = generative.Generative(drug_class)
# res, target_res, dock_target = generate.run()
# res_predict = {'drug_sequence' : res, 'drug_target' : target_res, 'drug_docking' : dock_target}

# export_dict_to_json(res_predict, f"{generate_random_code()}.json")
    
def DrugDeNovo(dc1="none", dc2="none", dc3="none", dc4="none", dc5="none", dc6="none", dc7="none", dir_target="."):
    drug_class = [dc1, dc2, dc3, dc4, dc5, dc6, dc7]
    drug_class = [drug_class_dict[dc] for dc in drug_class]
    print(drug_class)
    
    try:
      generate = generative.Generative(drug_class)
    except:
      generate = AIGDrugDeNovo.generative.Generative(drug_class)
    res, target_res, dock_target = generate.run()
    res_predict = {'drug_sequence' : res, 'drug_target' : target_res, 'drug_docking' : dock_target}
    filename_json = f"{dir_target}/{generate_random_code()}.json"
    export_dict_to_json(res_predict, filename_json)
    return filename_json
