var d3 = require('d3');
var $ = require('jquery');

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
