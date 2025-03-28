from xml.etree import ElementTree as ET
import xml.dom.minidom as minidom
from typing import List, Dict, Any, Union

def create_mitsuba_xml(
    output_file: str,
    shapes: List[Dict[str, Any]],
    materials: List[Dict[str, Union[str, float]]],
) -> str: 
    """
    Create a Mitsuba 3.6.2 compatible XML file.
    
    Parameters:
    -----------
    output_file : str
        Path to save the XML file
    shapes : list of dict
        List of shape dictionaries with properties like id, filename, material_id
    materials : list of dict
        List of material dictionaries with properties like id, type, thickness
    """
    # Create the root element
    root = ET.Element("scene")
    root.set("version", "2.1.0")

    # Add a comment for materials section
    materials_comment = ET.Comment(" Materials ")
    root.append(materials_comment)

    # Add materials
    for material in materials:
        bsdf = ET.SubElement(root, "bsdf")
        bsdf.set("type", "itu-radio-material")
        bsdf.set("id", material["id"])

        type_elem = ET.SubElement(bsdf, "string")
        type_elem.set("name", "type")
        type_elem.set("value", material["type"])

        thickness_elem = ET.SubElement(bsdf, "float")
        thickness_elem.set("name", "thickness")
        thickness_elem.set("value", str(material["thickness"]))

    # Add a comment for shapes section
    shapes_comment = ET.Comment(" Shapes ")
    root.append(shapes_comment)

    # Add shapes
    for shape in shapes:
        shape_elem = ET.SubElement(root, "shape")
        shape_elem.set("type", "ply")
        shape_elem.set("id", shape["id"])

        filename_elem = ET.SubElement(shape_elem, "string")
        filename_elem.set("name", "filename")
        filename_elem.set("value", shape["filename"])

        normals_elem = ET.SubElement(shape_elem, "boolean")
        normals_elem.set("name", "face_normals")
        normals_elem.set("value", "true")

        ref_elem = ET.SubElement(shape_elem, "ref")
        ref_elem.set("id", shape["material_id"])
        ref_elem.set("name", "bsdf")
        
        if "to_world" in shape:
            to_world_elem = ET.SubElement(shape_elem, "transform")
            to_world_elem.set("name", "to_world")
            to_world_string = " ".join([f"{x:.6f}" for x in shape["to_world"]])
            to_world_elem.set("value", to_world_string)

    # Convert to a string with pretty formatting
    rough_string = ET.tostring(root, encoding='utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    # Remove extra blank lines that minidom adds
    pretty_lines = [line for line in pretty_xml.split('\n') if line.strip()]
    clean_xml = '\n'.join(pretty_lines)

    # Write to file
    with open(output_file, 'w') as f:
        f.write(clean_xml)

    print(f"Created Mitsuba XML file: {output_file}")
    return clean_xml


def get_transformation_matrix(translation: List[float]=[0,0,0], scale: Union[float,List[float]]=1.0)->List[float]:
    """
    Get the transformation matrix for translation and scaling.
    
    Parameters:
    -----------
    translation : list of float, optional
        List of 3 floats for translation along x, y, z axes, by default [0,0,0]
    scale : float or list of float, optional
        Scaling factor for uniform scaling or list of 3 floats for non-uniform scaling, by default 1.0
    
    Returns:
    --------
    list of float
        List of 16 floats representing the transformation matrix
    """
    if isinstance(scale, float) or isinstance(scale, int):
        scale = [scale, scale, scale]
    return [
        scale[0], 0, 0, translation[0],
        0, scale[1], 0, translation[1],
        0, 0, scale[2], translation[2],
        0, 0, 0, 1
    ]