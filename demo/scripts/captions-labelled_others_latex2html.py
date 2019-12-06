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

def replace_caption_syntax(captions_matched_raw, cell_number, item):
    caption = captions_matched_raw[0]
    cell_content = nb['cells'][cell_number]['source']
    labels_matched_raw = [x for x in cell_content if re.search(labels_pattern_html, x)]
    #############################################
    # Here the label pattern can be modified either as LaTeX or HTML
    #############################################
    labels_tags_dic = get_labels_tags_dic(labels_pattern_html, nb)
    label = re.findall(labels_pattern_html, labels_matched_raw[0])
    #############################################
    html_syntax_items = (
        '<figcaption style="text-align:center;font-size:14px">',
        '<b>',
        marker+':'+labels_tags_dic[label[0]]+' ',
        '</b>',
        '<em>'+' ',
        caption,
        '</em>',
        '</figcaption>'
    )
    chunk_new = ''.join(html_syntax_items)
    chunk_original = r'\caption{'+caption+'}'
    item_modified = item.replace(chunk_original, chunk_new)
    return item_modified

def latex_to_html_captions(nb):
    for i in range(len(nb['cells'])):
        cell = nb['cells'][i]
        for j in range(len(cell['source'])):
            item = cell['source'][j]

            captions_matched_raw = re.findall(captions_pattern,item)
            if captions_matched_raw:
                cell_number = i
                item_modified = replace_caption_syntax(captions_matched_raw, cell_number, item)
                nb['cells'][i]['source'][j] = item_modified
    return nb

if __name__ == "__main__":
    import re
    import os
    import sys
    import json
    #############################################
    # Script input arguments (self-explanatory)
    #############################################
    nb_infile = sys.argv[1]
    nb_outfile = sys.argv[2]
    marker = sys.argv[3]
    #############################################
    # Settings
    #############################################
    strings_for_removal = ').,;:'
    labels_pattern_html = re.compile(r'<a id=(.+?)></a>')
    captions_pattern = re.compile(r'\\caption\{(.+?)\}$')
    #############################################
    # Load Jupyter Notebook as Dictionary
    #############################################
    try:
        print(os.path.isfile(nb_infile))
        with open(nb_infile,'r') as fp:
            nb = json.load(fp)
        fp.close()
#         print(nb)
        #############################################
        # Turn LaTeX label/ref syntax in Jupyter Noteboks into MD syntax
        #############################################
#         nb_modified = latex_to_html_captions(nb)
        
        print(nb)
        with open(nb_outfile,'w') as fp:
            json.dump(nb_modified, fp)
        fp.close()
        print("New Notebook successfully generated")
    except:
        print("Couldn't find Jupyter Notebook. Check your input path")
