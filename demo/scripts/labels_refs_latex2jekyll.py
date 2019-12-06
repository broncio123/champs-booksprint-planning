
def minimal_ordered_set(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def get_labels_tags_dic(labels_pattern, nb):
    COUNTER = 0
    labels_tags_all = []
    
    for i in range(len(nb['cells'])):
        cell = nb['cells'][i]
        for j in range(len(cell['source'])):
            item = cell['source'][j]
            labels_matched_raw = re.findall(labels_pattern, item)

            if labels_matched_raw:
                COUNTER +=1
                labels_tags = [s.translate({ord(i): None for i in strings_for_removal}) for s in labels_matched_raw]
                labels_tags = list(set(labels_tags))
                labels_tags_all += labels_tags
    #############################################
    # Test for non-repeated labels
    #############################################
    if COUNTER == len(labels_tags_all):
        print("Good! There aren't repeated labels in your file")
        labels_tags_dic = {labels_tags_all[n]:str(n+1) for n in range(COUNTER)}
        return labels_tags_dic
    else:
        print("There are repeated labels in your file. This will lead to errors")

def replace_ref_syntax(refs_matched_raw, item):
    refs_tags = [s.translate({ord(i): None for i in strings_for_removal}) for s in refs_matched_raw]
    refs_tags = list(set(refs_tags))
    #############################################
    # Replace Markdown reference syntax
    #############################################
    x = refs_tags[0]
    chunk_original = r'\ref{'+marker+':'+x+'}'
    x_number = labels_tags_dic[x]
    chunk_new = '['+marker+':'+x_number+'](#'+x+')'
    item_modified = item.replace(chunk_original, chunk_new)
    return item_modified
    
def replace_label_syntax(labels_matched_raw, item):
    labels_tags = [s.translate({ord(i): None for i in strings_for_removal}) for s in labels_matched_raw]
    labels_tags = list(set(labels_tags))
    #############################################
    # Replace Markdown reference syntax
    #############################################
    x = labels_tags[0]
    chunk_original = r'\label{'+marker+':'+x+'}'
    chunk_new = '\n<a id="'+x+'"></a>'
    item_modified = item.replace(chunk_original, chunk_new)
    return item_modified

def latex_to_jekyll_labels_refs(nb):
    #############################################
    # Search for LaTeX label/ref commands within cells' content (source)
    #############################################
    for i in range(len(nb['cells'])):
        cell = nb['cells'][i]
        for j in range(len(cell['source'])):
            item = cell['source'][j]
            refs_matched_raw = re.findall(refs_pattern,item)
            labels_matched_raw = re.findall(labels_pattern, item)
            #############################################
            # Modify cell content with new syntax accordingly
            #############################################
            if refs_matched_raw and not labels_matched_raw:
                item_modified = replace_ref_syntax(refs_matched_raw, item)
                nb['cells'][i]['source'][j] = item_modified

            elif labels_matched_raw and not refs_matched_raw:
                item_modified = replace_label_syntax(labels_matched_raw, item)
                nb['cells'][i]['source'][j] = item_modified

            else:
                pass
    #############################################
    # Modify cell content with new syntax accordingly
    #############################################
    return nb 

if __name__ == "__main__":
    import re
    import sys
    import json
    #############################################
    # Settings
    #############################################
    marker = 'fig'
    strings_for_removal = ').,;:'

    labels_pattern = re.compile(r'\\label\{'+marker+'(.+?)\}')
    refs_pattern = re.compile(r'\\ref\{'+marker+'(.+?)\}')
    #############################################
    # Script input arguments 
    #############################################
    nb_infile = sys.argv[1] 
    nb_outfile = sys.argv[2]
    #############################################
    # Load Jupyter Notebook as Dictionary
    #############################################
    try:
        with open(nb_infile,'r') as fp:
            nb = json.load(fp)
        fp.close()
        #############################################
        # Turn LaTeX label/ref syntax in Jupyter Noteboks into MD syntax
        #############################################
        labels_tags_dic = get_labels_tags_dic(labels_pattern, nb)
        nb_modified = latex_to_jekyll_labels_refs(nb)
        
        with open(nb_outfile,'w') as fp:
            json.dump(nb_modified, fp)
        fp.close()
        print("New Notebook successfully generated")
    except:
        print("Couldn't find Jupyter Notebook. Check your input path")
