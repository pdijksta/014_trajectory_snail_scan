import subprocess

def send_to_desy_elog(author, title, severity, text, elog, image=None):
    """
    Send information to a supplied electronic logbook.
    Author: Christopher Behrens (DESY)
    """

    # The DOOCS elog expects an XML string in a particular format. This string
    # is beeing generated in the following as an initial list of strings.
    succeded = True  # indicator for a completely successful job
    # list beginning
    elogXMLStringList = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<entry>']

    # author information
    elogXMLStringList.append('<author>')
    elogXMLStringList.append(author)
    elogXMLStringList.append('</author>')
    # title information
    elogXMLStringList.append('<title>')
    elogXMLStringList.append(title)
    elogXMLStringList.append('</title>')
    # severity information
    elogXMLStringList.append('<severity>')
    elogXMLStringList.append(severity)
    elogXMLStringList.append('</severity>')
    # text information
    elogXMLStringList.append('<text>')
    elogXMLStringList.append(text)
    elogXMLStringList.append('</text>')
    # image information
    if image:
        try:
            # encodedImage = base64.b64encode(image)
            elogXMLStringList.append('<image>')
            elogXMLStringList.append(image)
            elogXMLStringList.append('</image>')
        except:  # make elog entry anyway, but return error (succeded = False)
            succeded = False
    # list end
    elogXMLStringList.append('</entry>')
    # join list to the final string
    elogXMLString = '\n'.join(elogXMLStringList)
    # open printer process
    try:
        lpr = subprocess.Popen(['/usr/bin/lp', '-o', 'raw', '-d', elog],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # send printer job
        lpr.communicate(elogXMLString.encode('utf-8'))
    except:
        succeded = False
    return succeded
