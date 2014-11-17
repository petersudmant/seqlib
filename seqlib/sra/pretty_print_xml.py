import xml.dom.minidom

xml = xml.dom.minidom.parse("SraExperimentPackage.xml") # or xml.dom.minidom.parseString(xml_string)
pretty_xml_as_string = xml.toprettyxml()
print pretty_xml_as_string
