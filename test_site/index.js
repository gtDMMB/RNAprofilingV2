
var selected_node_id = -1;
var highlighted_nodes = [];
var highlighted_radial_diagrams = [];
var svgDoc;
var leafList;
var nodeShapes = {};
var indepNodeData = {};
var indepNodeHelixData = {};
var nodeIndexData = {};
var bpCoverageData = {};

// https://stackoverflow.com/questions/13405129/javascript-create-and-save-file/30832210#30832210
// Function to download data to a file
function download(data, filename, type) {
    var file = new Blob([data], {type: type});
    if (window.navigator.msSaveOrOpenBlob) // IE10+
        window.navigator.msSaveOrOpenBlob(file, filename);
    else { // Others
        var a = document.createElement("a"),
                url = URL.createObjectURL(file);
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(function() {
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);  
        }, 0); 
    }
}

function nodeClickEffect(id, child) {
    //this.classList.add("bar-highlight");
    //console.log(this);

    return function() {
        //console.log(id);
        //var id = this.getAttributeNS(null,"id");
        if (id === selected_node_id)
            return

        var diagram_file = "arc_diagram/arc_diagram_" + id + ".svg"

        var diagram_image = document.getElementById("arc-diagram")
        diagram_image.src = diagram_file

        //tree.findOne("#"+id).fill("green");
        child.fill("#AAA");

        highlighted_nodes.map(function(node) {
            node.fill("none");
        });

        highlighted_nodes = [child];

        highlighted_radial_diagrams.map(function(node) {
            node.classList.remove("radial-highlight");
        });

        if (leafList.includes(id)) {
            let radial_holder = document.getElementById("radial_" + id);
            radial_holder.classList.add("radial-highlight");
            highlighted_radial_diagrams = [radial_holder];

            let column = document.getElementById("leaf-image-column");
            let hasVerticalScrollbar = column.scrollHeight > column.clientHeight;
            if (hasVerticalScrollbar)
                radial_holder.scrollIntoView();
        }

        if (id in indepNodeData)
        {
            var htmlTable = document.getElementById("indep-table");
            while (htmlTable.rows.length > 0)
                htmlTable.deleteRow(-1);

            {
                let tableData = indepNodeData[id]
                for (let row = 0; row < tableData.length; row++)
                {
                    let tr = htmlTable.insertRow(-1);
                    for (let col = 0; col < tableData[row].length; col++) {
                        if (tableData.length == 2 && col == 0)
                            var tableElement = document.createElement("th");
                        else if (row == 0 && tableData.length > 2)
                            var tableElement = document.createElement("th");
                        else if (tableData[0][col] == "")
                            var tableElement = document.createElement("th");
                        else if (tableData[row][0] == "")
                            var tableElement = document.createElement("th");
                        else
                            var tableElement = document.createElement("td");
                        tableElement.innerHTML = tableData[row][col];

                        tr.appendChild(tableElement);
                    }
                }
            }

            document.getElementById("indep-table-div").classList.remove("hidden");
            //document.getElementById("indep-feat-div").classList.remove("hidden");
        }
        else
        {
            document.getElementById("indep-table-div").classList.add("hidden");
            //document.getElementById("indep-feat-div").classList.add("hidden");
        }

        selected_node_id = id

        document.getElementById("download-frequent").disabled = false;
        document.getElementById("download-all").disabled = false;

    }
}

function diagramClickEvent(id, treeNode) {

    return function() {
        if (id === selected_node_id)
            return

        var diagram_file = "arc_diagram/arc_diagram_" + id + ".svg"

        var diagram_image = document.getElementById("arc-diagram")
        diagram_image.src = diagram_file

        //tree.findOne("#"+id).fill("green");
        treeNode.fill("#AAA");

        highlighted_nodes.map(function(node) {
            node.fill("none");
        });

        highlighted_nodes = [treeNode];

        highlighted_radial_diagrams.map(function(node) {
            node.classList.remove("radial-highlight");
        });

        let radial_holder = document.getElementById("radial_" + id);
        radial_holder.classList.add("radial-highlight");
        highlighted_radial_diagrams = [radial_holder];

        selected_node_id = id

        document.getElementById("download-frequent").disabled = false;
        document.getElementById("download-all").disabled = false;
    }
}

function buildFeatureLegendTable() {

    let table = document.getElementById("feature-legend-table");
    let data = JSON.parse(legendJSON);
    let headerKeys = Object.keys(data[0]);

    let tr_header = table.insertRow(-1);

    for (let idx = 0; idx < headerKeys.length; idx++) {
        
        let headerElement = document.createElement("th");
        headerElement.innerHTML = headerKeys[idx];

        tr_header.appendChild(headerElement);
    }

    for (let idx = 0; idx < data.length; idx++) {
        let tr = table.insertRow(-1);

        for (let jdx = 0; jdx < headerKeys.length; jdx++) {
            let tableElement = document.createElement("td");
            tableElement.innerHTML = data[idx][headerKeys[jdx]]

            tr.appendChild(tableElement);
        }
    }
}

function buildRawLegendTable() {

    let table = document.getElementById("raw-legend-table");
    let data = JSON.parse(additionalLegendJSON);
    if (data.length == 0)
    {
        document.getElementById("raw-legend-div").classList.add("hidden");
        return;
    }

    let headerKeys = Object.keys(data[0]);

    let tr_header = table.insertRow(-1);

    for (let idx = 0; idx < headerKeys.length; idx++) {
        
        let headerElement = document.createElement("th");
        headerElement.innerHTML = headerKeys[idx];

        tr_header.appendChild(headerElement);
    }

    for (let idx = 0; idx < data.length; idx++) {
        let tr = table.insertRow(-1);

        for (let jdx = 0; jdx < headerKeys.length; jdx++) {
            let tableElement = document.createElement("td");
            tableElement.innerHTML = data[idx][headerKeys[jdx]]

            tr.appendChild(tableElement);
        }
    }
}

function buildRadialDiagramColumn() {
    let diagramType = document.getElementById("diagram-type");

    let data = JSON.parse(leafJSON);
    for (let idx = 0; idx < data.length; idx++) {
        
        let image_holder = document.createElement("div");
        let image = document.createElement("img");
        let id = data[idx]["id"]
        if (diagramType.checked) 
            image.src = "radial_diagram/radial_diagram_" + id + ".svg";
        else
            image.src = "arc_diagram/arc_diagram_" + id + ".svg";
        image.setAttribute("width", "100%");
        image_holder.setAttribute("id", "radial_" + id);
        image_holder.classList.add("radial-diagram");
        document.getElementById("leaf-image-column").appendChild(image_holder);
        image_holder.appendChild(image);

        image.addEventListener("click", diagramClickEvent(id, nodeShapes[id]));
    }

}

function onDiagramTypeChange() {
    let diagramType = document.getElementById("diagram-type");

    let data = JSON.parse(leafJSON);
    for (let idx = 0; idx < data.length; idx++) {
        
        let id = data[idx]["id"]
        let image_holder = document.getElementById("radial_" + id);
        let image = image_holder.children[0];

        if (diagramType.checked) 
            image.src = "radial_diagram/radial_diagram_" + id + ".svg";
        else
            image.src = "arc_diagram/arc_diagram_" + id + ".svg";
    }
}

function singleStructButtonPress() {
    if (selected_node_id === -1)
        return

    let index = nodeIndexData[selected_node_id]["mostCommonIndex"];
    let data = [sampleGTBOLTZ[index]];
    let dataString = data.join("\n");

    download(
        dataString, 
        sequenceName + "_node_" + selected_node_id + "_most_frequent_structure.gtboltz", 
        "gtboltz");
}

function allStructButtonPress() {
    if (selected_node_id === -1)
        return

    let indices = nodeIndexData[selected_node_id]["allIndices"];
    let data = indices.map(function(index) { return sampleGTBOLTZ[index]; });
    let dataString = data.join("\n");

    download(
        dataString, 
        sequenceName + "_node_" + selected_node_id + "_all_structures.gtboltz", 
        "gtboltz");
}

function indepFeatButtonPress() {
    if (selected_node_id === -1)
        return

    let data = indepNodeHelixData[selected_node_id];
    let dataString = data.join("\n");

    download(
        dataString, 
        sequenceName + "_node_" + selected_node_id + "_indep_region_helix_classes.csv", 
        "gtboltz");
}

document.addEventListener('DOMContentLoaded', function() {

    console.log("Test2");
    var treediv = document.getElementById("treediv");

    var draw = SVG().addTo(treediv);
    var tree = draw.group();
    tree.svg(treeTXT);
    var treeWidth = tree.width();
    //scale_factor = Math.min(treediv.clientWidth / treeWidth, 1.);
    scale_factor = treediv.clientWidth / treeWidth;
    //tree.scale(scale_factor,scale_factor, origin=[0,0]);
    tree.transform({scale:scale_factor, origin:[0,0]});
    //tree.translate(4 * scale_factor, 1080 * scale_factor);

    graph0 = document.getElementById("graph0");
    svgRoot = graph0.parentElement;

    draw.width("100%");
    draw.height(scale_factor * svgRoot.height.baseVal.value);
    //draw.height("300%");

    docHeight = svgRoot.viewBox.baseVal.height;

    console.log(draw);

    console.log("Loaded");
    console.log(tree);

    Array.prototype.map.call(
        tree.find(".node"),
        function(node) {
            let id = node.attr("id");
            if (!isNaN(id)){
                //console.log(node.children())
                node.children().map(function(child) {
                    if (child.type === "ellipse" || child.type == "polygon") {
                        nodeShapes[id] = child;
                        //console.log(id);
                        bbox = child.bbox();
                        //console.log(bbox);
                        //console.log(child);
                        //console.log(draw.height())
                        let rect = draw.rect(bbox.w,bbox.h);
                        rect.move(bbox.x, bbox.y );
                        rect.opacity(0);
                        rect.addTo(graph0);

                        rect.on("click", nodeClickEffect(id, child));
                        //child.click(mouseOverEffect(id, tree));
                    }
                });
                //node.on("mouseover", mouseOverEffect(id));
            }
        });
    console.log("Test");

    tree.data = tree.data;

    buildFeatureLegendTable();
    buildRawLegendTable();

    window.addEventListener("resize", function() {
        scale_factor = Math.min(treediv.clientWidth / treeWidth, 1.);
        tree.transform({scale:scale_factor, origin:[0,0]});
        draw.height(scale_factor * svgRoot.height.baseVal.value);
    });

    buildRadialDiagramColumn();

    let leafData = JSON.parse(leafJSON);
    leafList = leafData.map(function(row) { return row["id"]; });

    document.getElementById("diagram-type").addEventListener("change", onDiagramTypeChange);

    indepNodeData = JSON.parse(indepNodeJSON);
    //indepNodeHelixData = JSON.parse(indepNodeHelixJSON);
    nodeIndexData = JSON.parse(nodeSampleIndicesJSON);
    //bpCoverageData = JSON.parse(bpCoverageJSON);

    sampleGTBOLTZ = sampleGTBOLTZ.split("\n")
        .filter(function(elem) { return (elem !== ""); });

    singleStructButton = document.getElementById("download-frequent");
    allStructButton = document.getElementById("download-all");

    singleStructButton.disabled = true;
    allStructButton.disabled = true;

    singleStructButton.addEventListener("click", singleStructButtonPress);
    allStructButton.addEventListener("click", allStructButtonPress);

    document.getElementById("sequence-name").innerHTML = sequenceName;
    document.title = sequenceName + " Profile"

});
