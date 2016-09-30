/*
 The sortable grid is based on code from https://github.com/irusri/correlation-explorer
 ("An interactive way to explore 2d arrays")

 */
import * as d3 from 'd3';
import $ from 'jquery';

// CORREXPLORER.JS
var stringify = function (x) {
    if (typeof(x) === 'number' || x === undefined) {
        return String(x);
        // otherwise it won't work for:
        // NaN, Infinity, undefined
    } else {
        return JSON.stringify(x);
    }
};

// to update number next to zoom slider
function outputUpdate(zoomvalue) {
    $('#zoom_level').val(zoomvalue);
}

var nugenes = false;

window.onload = function(){

    ingest();


    d3.select("button#sort").on("click", function() {
        $('#loader').show();
        d3.select("#right_col svg").remove();
        document.getElementById("zoom").value = 1.0;
        document.getElementById("zoom_level").value = 1.0;
        if(nugenes) {
            readFormValues();
            ingest();
        }
        else {
            makeMatrix();
            main(mtrx, sampleList, geneList);
        }

    });

    $('.sorter').change( function() {
        $('#sort').addClass('button-primary');
    });

    $('#genes').change( function() {
        nugenes = true;
    })

};

var geneProps = [ //list of gene-specific data items, to go to gene name mouseover
    "Gene","Uniprot","Biotype","Chr","Cytoband","Druggable","Cancer_Gene_Census","Oncogene","Tumor_suppressor","Drugs"
];

var main = function(corr, label_col, label_row){ //console.log("corr", corr, label_col, label_row);

    var transition_time = 1500;

    var renderarea = d3.select('div#right_col');

    var tooltip = renderarea.select('div.tooltip')
        .style("opacity", 1e-6);

    var svg = renderarea.append('svg')
        .attr('width', 3000)
        .attr('height', 3000);

    var sort_process = $("select#sort_func").val();
    d3.select("select#sort_func").on("change", function() {
        sort_process = this.value;
        reorder_matrix(last_k, last_what);
    });


    var row = corr;
    var col = d3.transpose(corr);


    // converts a matrix into a sparse-like entries
    // maybe 'expensive' for large matrices, but helps keeping code clean
    var indexify = function(mat){
        var res = [];
        for(var i = 0; i < mat.length; i++){
            for(var j = 0; j < mat[0].length; j++){
                res.push({i:i, j:j, val:mat[i][j]});
            }
        }
        return res;
    };

    var corr_data = indexify(corr); //console.log("corr_data", corr_data);
    var order_col = d3.range(label_col.length + 1);
    var order_row = d3.range(label_row.length + 1);



    var scale = d3.scaleLinear()
        .domain([0, d3.min([50, d3.max([label_col.length, label_row.length, 4])])])
        .range([0, parseFloat($("input#zoom").val()) * 600]);

    var zoomfactor = parseFloat($("input#zoom").val());

    d3.select("input#zoom").on("change", function() {
        zoomfactor = parseFloat(this.value);

        scale = d3.scaleLinear()
            .domain([0, d3.min([50, d3.max([label_col.length, label_row.length, 4])])])
            .range([0, parseFloat(this.value) * 600]);

        refresh_order();
    });

    var div = d3.select("body").append("div")
        .attr("class", "tool_tip")
        .style("opacity", 0);

    var dot_rad = 0.25; // radius for samples without CNA

    var label_space = 125;
    // I will make it also a function of scale and max label length

    var matrix = svg.append('g')
        .attr('class','matrix')
        .attr('transform', 'translate(' + (label_space + 10) + ',' + (label_space + 10) + ')');

    var pixel = matrix.selectAll('g.pixel').data(corr_data);


    // as of now, data not changable, only sortable
    pixel =  pixel.enter()
        .append('g')
        .attr('class', function (d) { return label_col[d.j] +' '+ label_row[d.i] +' pixel';})
        .attr('transform', function(d){
            return 'translate('+ scale(order_col[d.j]) + ',' + scale(order_row[d.i]) +') scale('+zoomfactor+')';
        })
        .merge(pixel);

    // vertical guides for z-score
    pixel.append("rect")
        .attr('class', 'guide')
        .attr('y', scale(0.43))
        .attr('x', scale(0))
        .attr('height', scale(0.06))
        .attr('width', scale(0.9));
    pixel.append("rect")
        .attr('class', 'guide')
        .attr('y', scale(0.2))
        .attr('x', scale(0))
        .attr('height', scale(0.06))
        .attr('width', scale(0.9));
    pixel.append("rect")
        .attr('class', 'guide')
        .attr('y', scale(0.66))
        .attr('x', scale(0))
        .attr('height', scale(0.06))
        .attr('width', scale(0.9));

    // add g element for focal mutation, store a bunch of variables
    var focus = pixel.append("g")
        .attr('class', function(d) { //console.log("CNVtype"allData[d][entry].CNV_type);
            var self = allData[label_col[d.j]][label_row[d.i]]; //console.log(self);
            self._focal = false; // variables for tracking glyph type
            self._score = self.RNA_z_score === "-" ? 0 : parseFloat(self.RNA_z_score);
            // check a few conditions and create attributes
            if (self.CNV_type == 'focal') { self._focal = true; }
            if (self._score > 2) {
                self._score = 2;
            }
            if (self._score < -2) {
                self._score = -2;
            }

            if (self.CNV_log2_max_adj === "-") { self._cnv = "normal"; } else {
                var r = parseFloat(self.CNV_log2_max_adj);
                if (r < -1) { self._cnv = "deep"; }
                else if (-1 <= r && r < -0.235) { self._cnv = "del"; }
                else if (0.2 >= r && r >= -0.235) { self._cnv = "normal"; }
                else if (0.584 >= r && r > 0.2) { self._cnv = "amp"; }
                else if (r > 0.584) { self._cnv = "high"; }
            }
            if (self._focal === true) {
                if (self.Variant_Impact !== '-') {
                    return self.Variant_Impact + ' focal';
                } else {
                    return 'focal';
                }
            } else {
                return 'hidden';
            }
        });

    pixel.on("mouseover", function(d) {
        var self = allData[label_col[d.j]][label_row[d.i]];
        var txt = "<p class='header'>Sample: <b>" + label_col[d.j] + "</b>, Gene: <b>" + label_row[d.i] + "</b></p>";

        function prop2txt(ob) {
            for (var propertyName in ob) {
                if (ob[propertyName] !== "-" && propertyName.charAt(0) !== "_" && geneProps.indexOf(propertyName) == -1) {
                    txt += propertyName + ": " + (ob[propertyName] == "[object Object]" ? "&nbsp;&nbsp;" : ob[propertyName] + "<br/>");
                    if (typeof ob[propertyName] === 'object') { prop2txt(ob[propertyName]); }
                }
            }
        }
        prop2txt(self);

        div.transition()
            .duration(200)
            .style("opacity", .95);
        div.html(txt)
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 150) + "px");
    })
        .on("mouseout", function(d) {
            div.transition()
                .duration(500)
                .style("opacity", 0);
        });

    // draw focal cross
    focus.append('path')
        .attr('d', d3.symbol()
            .size(90)
            .type(d3.symbolCross))
        .attr('transform', function(d){
            var self = allData[label_col[d.j]][label_row[d.i]];
            var tx = 'translate('+ scale(0.45) + ',';
            if (self.RNA_z_score !== '-') {
                return tx + scale(0.45 - (self._score * 0.25)/2) + ')';
            } else {
                return tx + scale(0.45) +')';
            }
        })
        .attr('class', function(d) {
            var self = allData[label_col[d.j]][label_row[d.i]];
            return (self.Variant_Impact !== '-' ? self.Variant_Impact : null);

        });

    // draw dot, vertical position indicating RNA z-score
    pixel.append('circle')
        .attr('class', function(d) {
            var self = allData[label_col[d.j]][label_row[d.i]];
            return (self.Variant_Impact !== '-' ? self.Variant_Impact : null);
        })
        .attr('cx', scale(0.45))
        // y coordinate = RNA_z_score
        .attr('cy', function(d) {
            var self = allData[label_col[d.j]][label_row[d.i]];
            if (self.RNA_z_score !== '-') {
                return scale(0.45 - (self._score * 0.25)/2);
            } else {
                return scale(0.45);
            }
        })
        .attr('r', scale(dot_rad));

    // mark CNVs
    pixel.append("path")
        .attr("class", "cnv")
        .attr("d", function(d) {
            var self = allData[label_col[d.j]][label_row[d.i]];
            var cx = scale(0.45);
            var cy = (self.RNA_z_score !== '-' ? scale(0.45 - (self._score * 0.25)/2) : scale(0.45));
            var start = 270;
            var end = 90;
            if(self._cnv == "high") {  } // use values set already
            if(self._cnv == "amp") { start = 315; end = 45; }
            if(self._cnv == "normal") { end = 270 }
            if(self._cnv == "del") { start = 135; end = 225;  }
            if(self._cnv == "deep") { start = 90; end = 270;  }

            return describeArc(cx, cy, scale(dot_rad), start, end);

        });

    // mark breakpoint
    pixel.append('circle')
        .attr('class', function(d){
            var self = allData[label_col[d.j]][label_row[d.i]];
            return (self.SV !== '-' ? "Breakpoint" : "hidden");
        })
        .attr('cx', scale(0.45))
        // y coordinate = RNA_z_score
        .attr('cy', function(d) {
            var self = allData[label_col[d.j]][label_row[d.i]];
            if (self.RNA_z_score !== '-') {
                return scale(0.45 - (self._score * 0.25)/2);
            } else {
                return scale(0.45);
            }
        })
        .attr('r', scale(dot_rad/2));

    // column labels: samples
    var tick_col = svg.append('g')
        .attr('class','ticks')
        .attr('transform', 'translate(' + (label_space + 10) + ',' + (label_space) + ')')
        .selectAll('text.tick')
        .data(label_col);

    tick_col =  tick_col.enter()
        .append('text')
        .attr('class','tick')
        .style('text-anchor', 'start')
        .attr('transform', function(d, i){return 'translate(' + parseInt(scale(order_col[i] + 0.7)) + ', 0), rotate(270)';})
        .attr('font-size', scale(0.7))
        .text(function(d){ return d; })
        // .on('mouseover', function(d, i){tick_mouseover(d, i, col[i], label_row);})
        // .on('mouseout', function(d){mouseout(d);})
        .on('click', function(d, i){reorder_matrix(i, 'col');})
        .merge(tick_col);

    // row labels: genes
    var tick_row = svg.append('g')
        .attr('class','ticks')
        .attr('transform', 'translate(' + (label_space) + ',' + (label_space + 10) + ')')
        .selectAll('text.tick')
        .data(label_row);

    tick_row = tick_row.enter()
        .append('text')
        .attr('class','tick')
        .attr('y', function(d, i){return scale(order_row[i] + 0.7);})
        .style('text-anchor', 'end')
        .attr('font-size', scale(0.7))
        //.attr('title','Click to sort by this gene; hover longer to see sorting scores.')
        .text(function(d){ return d; })
        .on('mouseover', function(d){
            var self = allData[sampleList[0]][d];
            var txt = "<p class='header'>Gene: <b>" + d + "</b></p>";

            function prop2txt(ob) {
                for (var propertyName in ob) {
                    if (ob[propertyName] !== "-" && geneProps.indexOf(propertyName) !== -1) { //console.log(propertyName, geneProps.indexOf(propertyName));
                        txt += propertyName + ": " + (ob[propertyName] == "[object Object]" ? "&nbsp;&nbsp;" : ob[propertyName] + "<br/>");
                        if (typeof ob[propertyName] === 'object') { prop2txt(ob[propertyName]); }
                    }
                }
            }
            prop2txt(self);

            div.transition()
                .delay(1500)
                .duration(200)
                .style("opacity", .95);
            div.html(txt)
                .style("left", (d3.event.pageX + 10) + "px")
                .style("top", (d3.event.pageY - 15) + "px");


        })
        .on("mouseout", function(d) {
            div.transition()
                .duration(500)
                .style("opacity", 0);
        })
        .on('click', function(d, i){reorder_matrix(i, 'row');})
        .merge(tick_row);

    // sort grid and labels
    var refresh_order = function(){ //console.log("refreshing order");
        tick_col.transition()
            .duration(transition_time)
            .attr('font-size', scale(0.7))
            .attr('transform', function(d, i){ return 'rotate(270 ' + scale(order_col[i] + 0.7) + ',0)';})
            .attr('x', function(d, i){ return scale(order_col[i] + 0.7); });

        tick_row.transition()
            .duration(transition_time)
            .attr('font-size', scale(0.7))
            .attr('y', function(d, i){return scale(order_row[i] + 0.7);});

        pixel.transition()
            .duration(transition_time)
            .attr('transform', function(d){
                return 'translate('+ scale(order_col[d.j]) + ',' + scale(order_row[d.i]) +') scale('+zoomfactor+')';
            });
        $('#loader').hide();
    };

// sorting algorithms
    var last_k = 0;
    var last_what = 'col';
    var reorder_matrix = function(k, what){
        last_k = k;
        last_what = what;
        var order = [];
        var vec = [];
        var labels = [];
        var vecs = [];
        if(what === 'row'){  //yes, we are sorting counterpart
            vec = row[k];
            vecs = row;
            labels = label_col;  //test is if it ok
        } else if ( what === 'col' ) {
            vec = col[k];
            vecs = col;
            labels = label_row;
        }
        var indices = d3.range(vec.length);
        switch (sort_process) {
            case "value":
                indices = indices.sort(function(a,b){return vec[b] - vec[a];});
                break;
            case "abs_value":
                indices = indices.sort(function(a,b){return Math.abs(vec[b]) - Math.abs(vec[a]);});
                break;
            case "original":
                break;
            case "alphabetic":
                indices = indices.sort(function(a,b){return Number(labels[a] > labels[b]) - 0.5;});
                break;
            case "similarity":
                // Ugly, but sometimes we want to sort the coordinate we have clicked
                // Also, as of now with no normalization etc
                indices = d3.range(vecs.length);
                indices = indices.sort(function(a,b){
                    var s = 0;
                    for(var i = 0; i < vec.length; i++){
                        s += (vecs[b][i] - vecs[a][i]) * vec[i];
                    }
                    return s;
                });
                if ( what === 'col' ){
                    order_col = reverse_permutation(indices);
                } //not else if!
                if ( what === 'row' ) {
                    order_row = reverse_permutation(indices);
                }
                refresh_order();
                return undefined;
        }
        if ( what === 'row' ){
            order_col = reverse_permutation(indices);
        } //not else if!
        if ( what === 'col' ) {
            order_row = reverse_permutation(indices);
        }
        refresh_order();
    };

    var reverse_permutation = function(vec){
        var res = [];
        for(var i = 0; i < vec.length; i++){
            res[vec[i]] = i;
        }
        return res;
    };

};


// ONCOPRINT.JS
var displayGenes = {},
    allData = {},
    data = [],
    mtrx = [],
    sampleList = [
        "G15-01T-D",
        "G15-03T-D",
        "G15-05T-D",
        "G15-06T-D",
        "G15-07T-D",
        "G15-09T-D",
        "G15-10T-D-run2",
        "G15-12T-D",
        "G15-13T-D",
        "G15-14T-D",
        "G15-15T-D",
        "G15-16T-D",
        "G15-17T-D",
        "G15-21T-D",
        "G15-23T-D",
        "G15-23T2-D",
        "G15-24T-D",
        "G15-25T-D",
        "G15-26T-D",
        "G15-27T-D",
        "G15-29T-D",
        "G15-31T1-D",
        "G15-31T2-D",
        "G15-32T-D",
        "G15-33T-D"
    ];
var geneList = [
    "ATR",
    "ATRX",
    "BRAF",
    "BRCA1",
    "BRCA2",
    "CDK4",
    "CDK6",
    "CDKN2A",
    "CDKN2B",
    "DNMT3A",
    "EGFR",
    "ERBB3",
    "ERG",
    "FGFR3",
    "FLT3",
    "HGF",
    "IDH1",
    "IDO2",
    "KDR",
    "KIT",
    "MDM2",
    "MET",
    "MGMT",
    "NF1",
    "NFKB1",
    "PDGFRA",
    "PHLPP2",
    "PIK3CA",
    "PIK3CG",
    "PIK3R1",
    "PSMD2",
    "PTEN",
    "RB1",
    "ROS1",
    "SMO",
    "STAG2",
    "SUFU",
    "TERT",
    "TP53"
];
var impact_mod, impact_hi, sv, cnv, rna; // variables for weighting values

export function ingest() {
    var samp = 0; // variable for keeping track of samples; it's incremented after data is loaded

    for (var j = 0; j < sampleList.length; j++) {
        allData[sampleList[j]] = {};

        //import txt file
        data[j] = d3.tsv("data/v2/" + sampleList[j] + ".txt",

            function(error, sampleData) {
                if (error) return console.warn(error);
                for (var n = sampleData.length - 1; n >= 0; n--) {
                    //  for each sample, look for genes from the geneList
                    //      for the next sample, do the same
                    //      if the GeneName is already present, add sample name to the list of values
                    //          e.g., {"GeneName1":["S1", "S2"]}

                    if (geneList.indexOf(sampleData[n].Gene) !== -1) {
                        //  add them to an object CCG ({"GeneName1":["S1"]})
                        if (!displayGenes[sampleData[n].Gene]) { displayGenes[sampleData[n].Gene] = []; }
                        displayGenes[sampleData[n].Gene].push(sampleList[samp]);

                        //add per-gene data\ to allData object
                        allData[sampleList[samp]][sampleData[n].Gene] = sampleData[n];
                    }
                }

                if (samp == sampleList.length - 1) {
                    console.log("displayGenes", displayGenes);
                    console.log("allData", allData);

                    makeMatrix();

                    $('#loader').hide();
                    main(mtrx, sampleList, geneList);

                }
                samp++;
            });
    }
}

function readFormValues() {
    impact_mod = parseFloat($("input#impact_mod").val());
    impact_hi = parseFloat($("input#impact_hi").val());
    sv = parseFloat($("input#sv").val());
    cnv = parseFloat($("input#cnv").val());
    rna = parseFloat($("input#rna").val());
    var gtemp = $("#genes").val().split(",");
    for (var i = gtemp.length - 1; i >= 0; i--) {
        gtemp[i] = gtemp[i].trim();
    }
    geneList = gtemp;
    if(nugenes) { console.log("new geneList", geneList); }
    console.log("weights", impact_mod, impact_hi, sv, cnv, rna);
}

function makeMatrix() {
    //reset
    mtrx = [];
    // define weighting inputs
    readFormValues();
    for (var g = 0; g < geneList.length; g++) {
        var row = [];
        for (var s = 0; s < sampleList.length; s++) {
            var thisGene = allData[sampleList[s]][geneList[g]];
            var val = 0;
            //console.log("thisGene", thisGene);
            if (thisGene.Variant_Impact === "MODERATE") { val += impact_mod; }
            if (thisGene.Variant_Impact === "HIGH") { val += impact_hi; }
            if (thisGene.SV !== "-") { val += sv; }
            if (thisGene.CNV_log2_max_adj !== "-") {
                val += parseFloat(thisGene.CNV_log2_max_adj) * cnv;
            }
            if (thisGene.RNA_z_score !== "-") {
                val += parseFloat(thisGene.RNA_z_score) * rna;
            }

            row.push(val);
        }
        mtrx.push(row);
    }
    $('#sort').removeClass('button-primary');
}

//functions for drawing the pie-shaped annotations for CNVs
function polarToCartesian(centerX, centerY, radius, angleInDegrees) {
    var angleInRadians = (angleInDegrees - 90) * Math.PI / 180.0;

    return {
        x: centerX + (radius * Math.cos(angleInRadians)),
        y: centerY + (radius * Math.sin(angleInRadians))
    };
}

function describeArc(x, y, radius, startAngle, endAngle) {

    var start = polarToCartesian(x, y, radius, endAngle);
    var end = polarToCartesian(x, y, radius, startAngle);
    var center = x + "," + y;

    var arcSweep = endAngle - startAngle <= 180 ? "0" : "1";


    var d = [
        "M", start.x, start.y,
        "A", radius, radius, 0, arcSweep, 0, end.x, end.y,
        "L", center, "Z"
    ].join(" ");

    return d;
}
// how to use: d3.select("#something").append("path").attr("d", describeArc(20, 40, 10, 0, 180));
