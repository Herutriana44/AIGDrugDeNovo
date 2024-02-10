import streamlit as st
import json
import generative

# Definisikan nama file json
filename = "drug_class_mapping.json"

# Buka file json
with open(filename, "r") as file:
    data = json.load(file)

# Buat form input
with st.form(key="drug_class_form"):
    drug_class1 = st.selectbox("Drug Class 1", list(data.keys()))
    drug_class2 = st.selectbox("Drug Class 2", list(data.keys()))
    drug_class3 = st.selectbox("Drug Class 3", list(data.keys()))
    drug_class4 = st.selectbox("Drug Class 4", list(data.keys()))
    drug_class5 = st.selectbox("Drug Class 5", list(data.keys()))
    drug_class6 = st.selectbox("Drug Class 6", list(data.keys()))
    drug_class7 = st.selectbox("Drug Class 7", list(data.keys()))

    submitted = st.form_submit_button("Submit")

# Jika tombol submit diklik
if submitted:
    drug_class = [drug_class1, drug_class2, drug_class3, drug_class4, drug_class5, drug_class6, drug_class7]
    drug_class = [data[dc] for dc in drug_class]
    # print(drug_class)
    generate = generative.Generative(drug_class)
    res, target_res, dock_target = generate.run()
    # print(res, target_res)
    st.caption("Drug Sequence")
    st.code(res, language="python")
    st.caption("Drug Target")
    st.code(target_res, language="python")
    st.caption("Drug Docking With Target")
    st.code(dock_target, language="python")
    # print(drug_class)

# generative = generative(drug_class)