import xml.etree.ElementTree as ET

def parse_sra_xml(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()

    biosamples = []
    for biosample in root.findall(".//BioSample"):
        biosample_id = biosample.find(".//Id[@db='BioSample']")
        sample_name = biosample.find(".//Id[@db_label='Sample name']")
        
        if biosample_id is not None and sample_name is not None:
            
            if 'Kronos' not in sample_name.text:
                biosamples.append((biosample_id.text, 'Kronos' + sample_name.text))
            else:
                biosamples.append((biosample_id.text, sample_name.text))
        else:
            print(f'{biosample_id} and {sample_name} not properly paired')

    return biosamples


file_path = 'C:\\Users\\skyun\\Downloads\\biosample_result.xml'  # Replace with the path to your XML file
biosamples = parse_sra_xml(file_path)

with open('C:\\Users\\skyun\\Downloads\\biosample2Kronos.list', 'w') as o:
    for item in biosamples:
        o.write('\t'.join(item) + '\n')
