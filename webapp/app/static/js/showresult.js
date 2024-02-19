/* Copyright (C) 2015,2016 Mohammad Alanjary
   University of Tuebingen
   Interfaculty Institute of Microbiology and Infection Medicine
   Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
   Funding by the German Centre for Infection Research (DZIF)

   This file is part of ARTS
   ARTS is free software. you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version

   License: You should have received a copy of the GNU General Public License v3 with ARTS
   A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.
*/

$('body').scrollspy({ target: '#spy', offset:100});

var summTable = 0;
var krTable = 0;
var dupTable = 0;
var bgcTable = 0;
var corecircle = 0;
var summchart = 0;
var statustimer = 0;
var genetreetimer = 0;
var finished = 0;
var all_regions = 0;

//Parses URI parameters
var getUrlParameter = function getUrlParameter(sParam) {
    var sPageURL = window.location.search.substring(1),
        sURLVariables = sPageURL.split('&'),
        sParameterName,
        i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return sParameterName[1] === undefined ? true : decodeURIComponent(sParameterName[1]);
        }
    }
};

//function rolloverpie(slice){
//    console.log("rollover",slice);
//}
//function rolloutpie(slice){
//    console.log("rollout",slice);
//}
var sortgeneselector = function() {
    var optarr = $("#genetreeselector option");
    optarr.sort(function(a,b) {
        if (a.text.toLowerCase() > b.text.toLowerCase()) return 1;
        if (a.text.toLowerCase() < b.text.toLowerCase()) return -1;
        return 0; })
    $("#genetreeselector").empty().append(optarr);
}

var d3pieobj = {"header":{"title":{"text":"Core Functions","fontSize":24,"font":"arial"},
                "subtitle":{"text":"","color":"#999999","fontSize":14,"font":"arial"},"titleSubtitlePadding":9},
                "size":{"canvasHeight":250,"canvasWidth":250,"pieInnerRadius":"50%","pieOuterRadius":"90%"},
                "data":{"sortOrder":"value-asc","content":[]},"labels":{"outer":{"format":"none","pieDistance":32},"inner":{"hideWhenLessThanPercentage":3},"mainLabel":{"fontSize":11},
                "percentage":{"color":"#ffffff","decimalPlaces":0},"value":{"color":"#adadad","fontSize":11},"lines":{"enabled":true},"truncation":{"enabled":true}},
//                "tooltips":{"enabled":true,"type":"caption","styles":{"fontSize":12}},
                "effects":{"pullOutSegmentOnClick":{"effect":"linear"}},
                "callbacks":{"onload":function(){
                    //remove oncick popout
                    var arcs = d3.selectAll(".cf_arc,.cf_labelGroup-inner,.cf_labelGroup-outer");
                    arcs.on("click",function(x){
                        $('#corecirclabel').html(x.label + " (" + x.value + ")");
                    });
                    arcs.on("mouseover",function(x){
                        $('#corecirclabel').html(x.label + " (" + x.value + ")");
                    });
                    arcs.on("mouseout",function(x){
                        $('#corecirclabel').html("");
                    });
                }},
                "misc":{"cssPrefix":"cf_","gradient":{"enabled":true,"percentage":100}}};


function updatestatus(data){
    if (data.orgname) {
        $('#jobtitle').html(data.jobtitle);
        $('#corecount').html(data.coretotal);
        $('#cdscount').html(data.cdscount);
        $('#bgccount').html(data.bgccount);
        $('#dupcount').html(data.dupcount);
        $('#proxcount').html(data.proxcount);
        $('#phylcount').html(data.phylcount);
        $('#twocount').html(data.twocount);
        $('#threecount').html(data.threecount);
        $('#knownhitscount').html(data.krhits);
        }
    if (!corecircle && data.coretotal > 0 && data.cdscount > 0) {
//        corecircle = $("#corecircle").circliful({
//        animationStep: 1,
//        foregroundBorderWidth: 5,
//        backgroundBorderWidth: 15,
//        percent: parseInt(100*data.coretotal/data.cdscount),
//        foregroundColor: "#AA3838",
//        fontColor: "#AA3838",
//        textColor: "#AA3838",
//        halfCircle: 1,
//        text: "Core Genes"
//        });
        resp = $.getJSON("funcstats",function(funcdata){

            var colors = d3.scale.category20().range();
            var i=0;
            var content=[];
//            var corekeycont = "<ul class='list-inline corecirckey'>"
            $.each(funcdata,function(key,value){
                caption = key;
                if (caption.length > 20) caption = caption.substring(0,20)+"...";
                if (key!="data") content.push({"label":key,"value":value,"color":colors[i],"caption":caption+" "+value});
//                corekeycont+="<li><div class='smallbox' style='background-color:" + colors[i] + "''></div>" + key + "</li>";
                i++;
            });
//            corekeycont+="</ul>"
//            $("#corekey").html(corekeycont);
            if (content.length >= 1){
                d3pieobj.data.content = content;
                d3pieobj.header.subtitle.text = parseInt(100*data.coretotal/data.cdscount)+"% of total genes";
                corecircle = new d3pie("corecircle",d3pieobj);
                }
            });
     }

    pwidth = data.pwidth || "5%";
    ptitle = data.ptitle || "";

    if (data.state) {
        $('#statetxt').html(data.state);
        if (data.state.toLowerCase() == "error"){
            $('#jobprogbar').addClass("progress-bar-danger");
            $('#jobprogbar').removeClass("progress-bar-striped");
            $('#jobprogbar').removeClass("active");
            ptitle = "Error at step: "+data.step+", see log"
            pwidth = "100%"
        }
        if (data.state.toLowerCase() == "done"){
            $('#jobprogbar').addClass("progress-bar-success");
            $('#jobprogbar').removeClass("progress-bar-striped");
            $('#jobprogbar').removeClass("active");
            pwidth = "100%"
        }
    }
    $('#jobprogbar').width(pwidth);
    $('#steptxt').html(ptitle);

    var date = new Date();
    $.get("log?_="+date.getTime(),function (data){
        $("#logtxt").html(data.replace(/\n/g, "<br />"));
    });

    if (data.end || data.state == "Done" || finished){
        finished = 1;
        clearInterval(statustimer);
        explinks = '<a class="btn btn-primary" style="width:400px" href="xlfile">Export Tables</a><br><br>';
        explinks += '<a class="btn btn-primary" style="width:400px" href="export/alltrees.zip">Download All Trees</a><br><br>';
        explinks += '<a class="btn btn-primary" style="width:400px" href="export/aligned_core_genes.zip">Download Aligned Core Gene Fastas</a><br><br>';
        explinks += '<a class="btn btn-primary" style="width:400px" href="/archive/'+$("#jobid").val()+'.zip">Zip all results (Viewable in webapp)</a>';
        $('#exportlinks').html(explinks);
        $('#sptreediv').html("<img src='trees/speciesmlst' id='sptreeimg' width='800'>");

        summTable.ajax.reload();
        krTable.ajax.reload();
        dupTable.ajax.reload();
        bgcTable.ajax.reload();

        setTimeout(function(){
            $('[data-toggle="tooltip"]').tooltip()
            },2000);
    }
    if (summTable.rows().count() == 0){
        summTable.ajax.reload();
    }
    if (krTable.rows().count() == 0){
        krTable.ajax.reload();
    }
    if (dupTable.rows().count() == 0){
        dupTable.ajax.reload();
    }
    if (bgcTable.rows().count() == 0){
        bgcTable.ajax.reload();
    }
//    if (!all_regions) {
//        var http = new XMLHttpRequest();
//        http.open('HEAD', "antismash/all_regions.js");
//        http.send();
//        if (http.status != 404) $.getScript("antismash/all_regions.js");
//    }
//    console.log("update status");
}

//$(window).focus(function() {
//    if (!statustimer && !finished){
//        location.reload();
//    }
//    if (!statustimer && !finished){
//        $.getJSON('status',updatestatus);
//        clearInterval(statustimer);
//        statustimer = setInterval(function() {
//        $.getJSON('status',updatestatus);
//        }, 5000);
//    }
//});
//
//$(window).blur(function() {
//    clearInterval(statustimer);
//    statustimer = 0;
//});

//DATATABLES
$(document).ready(function(){
    krTable = $('#krTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "krtab"
            },
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        columns: [
            {"data":0},
            {"data":1},
            {"data":2},
            {"data":3},
            {"data":4},
            {"data":5, "visible":false,"defaultContent":""},
        ],
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
        autoWidth: false
    });
    dupTable = $('#dupTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "dupmatrix"
            },
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        columns: [
            {
                "className":      'details-control',
                "orderable":      false,
                "data":           null,
                "defaultContent": '',
                "width":"25px"
            },
            {"data":0},
            {"data":7,"defaultContent":"n/a","visible":true},
            {"data":1},
            {"data":2},
            {"data":3},
            {"data":4},
            {"data":5},
            {"data":6,"visible":false,"defaultContent":""}
        ],
        createdRow: function ( row, data, index ) {
            $('td',row).eq(5).html((data[3]*100).toFixed(1)+"%")
            $('td',row).eq(6).html((data[4]*100).toFixed(1)+"%")
            $('td',row).eq(7).html((data[5]*100).toFixed(1)+"%")
        },
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
        autoWidth: false,
        order: [4,'desc']
    });
    bgcTable = $('#bgcTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "bgctable"
            },
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        columns: [
            {
                "className":      'details-control',
                "orderable":      false,
                "data":           null,
                "defaultContent": '',
                "width":"25px"
            },
            {"data":0},
            {"data":1},
            {"data":2},
            {"data":3},
            {"data":4},
            {"data":5}
        ],
        createdRow: function ( row, data, index ) {
            $(row).addClass("bgcrow");
            bgcID = data[0].replace("cluster-","r").replace("_","c");
            $('td',row).eq(1).html("<a href='antismash/index.html#"+bgcID+"' onclick='openaslink(this)' target='bgcwindow' class='aslink' id='"+bgcID+"-link' title='view in antismash' data-toggle='tooltip' data-placement='right'>"+data[0].replace("cluster-","r").replace("_","c")+"</a>");
            if (data[5].length == 0){
                $('td',row).eq(0).removeClass();
                $('td',row).eq(0).html("<p align='center'><i>n/a</i></p>");
            }
            else{
                $('td',row).eq(0).attr("title","Expand/collapse row");
            }
            if (!all_regions){
             all_regions=1;
             $.ajax({
                cache: false,
                url: "antismash/regions.js",
                success: null,
                error:function (xhr, ajaxOptions, thrownError){
                    if(xhr.status==404) {
                        console.log(thrownError);
                    }
                }
                });
            }
        },
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: -1,
        autoWidth: false,
        order: [5,'desc']
    });
    summTable = $('#summaryTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "summarytab"
            },
        columns: [
            {
                "className":      'details-control',
                "orderable":      false,
                "data":           null,
                "defaultContent": '',
                "width":"25px"
            },
            {"data":"coregene"},
            {"data":"description"},
            {"data":"func"},
            {"data":"TC","visible":false,"defaultContent":"n/a"},
            {"data":"dNdS","visible":false,"defaultContent":"n/a"},
            {"data":"SC","visible":false,"defaultContent":"n/a"},
            {"data":"Ubiq","visible":true,"defaultContent":"n/a"},
            {"data":"duplicate", "width":"10px", "defaultContent":"n/a"},
            {"data":"proximity", "width":"10px","defaultContent":"n/a"},
            {"data":"phylogeny", "width":"10px","defaultContent":"n/a"},
            {"data":"known_hit", "width":"10px","defaultContent":"n/a"},
        ],
        order: [[11, 'desc'],[9, 'desc']],
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        createdRow: function ( row, data, index ) {
                if (data.duplicate == "Yes"){
                        $('td',row).eq(5).html("<span class='glyphicon glyphicon-ok'></span>");
                        $('td',row).eq(5).addClass("summtab_yes");
                    }
                else if (data.duplicate == "No"){
                        $('td',row).eq(5).html("<span class='glyphicon glyphicon-remove'></span>");
                    }
                if (data.proximity == "Yes"){
                        $('td',row).eq(6).html("<span class='glyphicon glyphicon-ok'></span>");
                        $('td',row).eq(6).addClass("summtab_yes");
                        $("#asresultslink").html("<a href='antismash/index.html' class='btn btn-primary' target='bgcwindow'>View BGC overview</a> <a href='#' class='btn btn-primary' onclick='toggleallbgcrow(event)'>Expand/Collapse all</a>");
                    }
                else if (data.proximity == "No"){
                        $('td',row).eq(6).html("<span class='glyphicon glyphicon-remove'></span>");
                    }
                if (data.known_hit == "Yes"){
                        $('td',row).eq(8).html("<span class='glyphicon glyphicon-ok'></span>");
                        $('td',row).eq(8).addClass("summtab_yes");
                    }
                else if (data.known_hit == "No"){
                        $('td',row).eq(8).html("<span class='glyphicon glyphicon-remove'></span>");
                    }
                if (data.phylogeny == "Yes"){
                        $('td',row).eq(7).html("<span class='glyphicon glyphicon-ok'></span>");
                        $('td',row).eq(7).addClass("summtab_yes");
                        if ($("#genetreeselector option[value='"+data.coregene+"']").length == 0){
                            labeltext = data.description;
                            if (labeltext.length >= 30){
                                labeltext = labeltext.substring(0,29)+"...";
                            }
                            $("#genetreeselector").append($("<option></option>").attr("value",data.coregene).text(labeltext));
                            $("#gntreediv").html("<h5><span class='glyphicon glyphicon-arrow-up'></span> Select a gene tree above to view tree</h5>");
                        }
                        clearTimeout(genetreetimer);
                        genetreetimer = setTimeout(sortgeneselector,3);
                    }
                else if (data.phylogeny == "No"){
                        $('td',row).eq(7).html("<span class='glyphicon glyphicon-remove'></span>");
                    }

            },
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
        autoWidth: false
    });
    // Update page every 5 sec
    $.getJSON('status',updatestatus);
    clearInterval(statustimer);
    statustimer = setInterval(function() {
        $.getJSON('status',updatestatus)
    }, 5000);
});

function update_size(){
    $('#summaryTable').css({ width: $('#summaryTable').parent().width() });
    if (summTable && summTable.fnAdjustColumnSizing) {
        summTable.fnAdjustColumnSizing();
    }
    $('#krTable').css({ width: $('#krTable').parent().width() });
    if (krTable && krTable.fnAdjustColumnSizing) {
        krTable.fnAdjustColumnSizing();
    }
    $('#dupTable').css({ width: $('#dupTable').parent().width() });
    if (dupTable && dupTable.fnAdjustColumnSizing) {
        dupTable.fnAdjustColumnSizing();
    }
    $('#bgcTable').css({ width: $('#bgcTable').parent().width() });
    if (bgcTable && bgcTable.fnAdjustColumnSizing) {
        bgcTable.fnAdjustColumnSizing();
    }
  }

//$(window).resize(function() {
//    clearTimeout(window.refresh_size);
//    window.refresh_size = setTimeout(function() { update_size(); }, 250);
//  });
//$(window).resize(update_size());

/* Formatting function for row details */
function formatbgcrow ( d, clust ) {
    if (clust) {
        childstring = "<div class='bgcdiv' id='"+d[0].replace("cluster-","r").replace("_","c")+"-svg'></div><br>";
        childstring += '<div class="legend"><div><div class="legend-field legend-type-biosynthetic"></div><div class="legend-label">Biosynthetic</div></div><div><div class="legend-field legend-type-transport"></div>';
        childstring += '<div class="legend-label">Transport</div></div><div><div class="legend-field legend-type-regulatory"></div><div class="legend-label">Regulatory</div></div><div>';
        childstring += '<div class="legend-field legend-type-other"></div><div class="legend-label">Other</div><div class="legend-field legend-type-DUF"></div><div class="legend-label">DUF</div>';
        childstring += '<div class="legend-field legend-type-Core"></div><div class="legend-label">Core</div><div class="legend-field legend-type-ResModel"></div><div class="legend-label">ResModel</div><div class="legend-field legend-type-CoreRes"></div><div class="legend-label">CoreRes</div></div></div><br><br>';
        }
    else childstring = "<div style='text-align:center'>No cluster annotation found</div><br>";
    childstring += '<table cellpadding="5" class="table bgcChild">';
    childstring += '<tr><th>Sequence id</th><th>Location (start-end)</th><th>Type</th><th>Gene</th><th>Description</th><th>Function</th></tr>';
    for (i = 0; i < d[6].length; i++){
//        clust..push({start:parseInt(d[5][i][2]),strand:1,end:parseInt(d[5][i][3]),locus_tag:d[5][i][1],type:d[5][i][4],description:d[5][i][5]});
        if (clust) {
            for (j = 0; j <clust.orfs.length; j++){
                //Test for orf positions that are not exact. Crude check for overlap > 50%
                start = parseInt(d[6][i][2]);
                end = parseInt(d[6][i][3]);
                orflen = clust.orfs[j].end - clust.orfs[j].start;
                overlap = Math.max(0, Math.min(clust.orfs[j].end, end) - Math.max(clust.orfs[j].start, start));
                if ( overlap/orflen > 0.5 ) {
                    clust.orfs[j].locus_tag += ":"+d[6][i][1];
                    clust.orfs[j].description += " [ "+d[6][i][5]+" ]";
                    if ((clust.orfs[j].type == "ResModel" && d[6][i][4] == "Core") || (clust.orfs[j].type == "Core" && d[6][i][4] == "ResModel")) {
                        clust.orfs[j].type = "CoreRes";
                        d[6][i][4] = "CoreRes";
                    }
                    else{
                        clust.orfs[j].type = d[6][i][4];
                    }
                    break;
                }
            }
        }
        childstring += '<tr><td>'+d[6][i][0]+'</td><td>'+d[6][i][2]+' - '+d[6][i][3]+'</td><td class="hittype_'+d[6][i][4]+'">'+d[6][i][4]+'</td><td>'+d[6][i][1]+'</td><td>'+d[6][i][5]+'</td><td>'+d[6][i][6]+'</td></tr>';
    }
    childstring += "</table>"
    return [childstring,clust];
}

function formatchildrow ( d ) {
    // `d` is the original data object for the row
    summary = summTable.ajax.json()
    childstring = '<table cellpadding="5" class="table summtabChild">';
    childstring += '<tr><th>Sequence id</th><th>Source</th><th>Location (start-end)</th><th>Bit Score</th><th>BGC</th><th>Type</th><th>Phylogeny</th></tr>';
    for (i = 0; i < d.hits.length; i++){
        bgc = d.proxhits[d.hits[i]] || ["","","","",""];
        phyl = "";
        if (("phylhits" in d) && d.phylhits[d.hits[i]]) phyl = "<a href='#phylogeny' onclick='changegnselector("+'"'+d.coregene+'"'+")'>View tree</a>";

        seq = summary["seqs"][d.hits[i]];
        childstring += "<tr><td class='seqid'><span data-toggle='tooltip' data-placement='right' title='"+seq[6].replace(/'/g,"")+"'>"+d.hits[i]+'</span></td><td>'+seq[2]+"</td><td>"+ seq[3]+"-"+seq[4]+"</td><td>"+seq[7];
        if (bgc[0].length) {
            bgcId = 'r'+bgc[0].replace('_','c')
            childstring += "</td><td><a href='#' onclick=\"findbgc('"+bgcId+"');\">"+bgc[0]+"</a>";
            }
        else {
            childstring += "</td><td>"+bgc[0];
            }
        childstring += "</td><td>"+bgc[1]+"</td><td>"+phyl+"</td></tr>";
    }
    childstring += "</table>";
    return childstring;
}

    // Add event listener for opening and closing details
$('#summaryTable tbody').on('click', 'td.details-control', function () {
        var tr = $(this).closest('tr');
        var row = summTable.row( tr );
        if ( row.child.isShown() ) {
            // This row is already open - close it
            row.child.hide();
            tr.removeClass('shown');
        }
        else {
            // Open this row
            row.child( formatchildrow(row.data()) ).show();
            tr.addClass('shown');
            $('[data-toggle="tooltip"]').tooltip();
        }
    } );

function formatduprow( hitlist ){
    childstring = "<ul>";
    for (i = 0; i < hitlist.length; i++){
        childstring += "<li>"+hitlist[i]+"</li>";
    }
    childstring += "</ul>";
    return childstring;
}

$('#dupTable tbody').on('click', 'td.details-control', function () {
        var tr = $(this).closest('tr');
        var row = dupTable.row( tr );
        if ( row.child.isShown() ) {
            // This row is already open - close it
            row.child.hide();
            tr.removeClass('shown');
        }
        else {
            // Open this row
            row.child( formatduprow(dupTable.ajax.json()["hits"][row.data()[0]]) ).show();
            tr.addClass('shown');
            $('[data-toggle="tooltip"]').tooltip();
        }
    } );

function openbgctr(tr){
    var row = bgcTable.row( tr );
    d = row.data();
    clust = all_regions[d[0].replace("cluster-","r").replace("_","c")];
    rslt = formatbgcrow(d,clust);
    row.child( rslt[0] ).show();
    tr.addClass('shown');
    if (clust) svgene.drawClusters(d[0].replace("cluster-","r").replace("_","c")+"-svg", [rslt[1]], 20, 800);
}
function closebgctr(tr){
    var row = bgcTable.row( tr );
    row.child.hide();
    tr.removeClass('shown');
}
function togglebgctr(tr){
    var row = bgcTable.row( tr );
    if ( row.child.isShown()) {
        closebgctr(tr)
    }
    else {
        // Open this row
        openbgctr(tr)
    }
}

$('#bgcTable tbody').on('click', 'td.details-control', function () {
        var tr = $(this).closest('tr',false);
        togglebgctr(tr);
    } );

//function exporttablesbtn(){
//    $("a.dt-button").each(function clickexport(){
//        $(this).click();
//    });
//}

$('a.zoombtn').click(function zoomimg(e){
    e.preventDefault();
    imgid = $(this).data("imgid");
    zoom = $(this).data("zoom");
    imgw = $("#"+imgid).width();
    if (zoom == "out" && imgw > 150){
        $("#"+imgid).width(imgw - 40);
    }
    if (zoom == "in" && imgw < 1500){
        $("#"+imgid).width(imgw + 40);
    }
});

function openbgcrow(x,y){
    y = typeof y !== 'undefined' ? y : 0;
    $("#bgcTable > tbody > tr.bgcrow").each(function(i,tr){
        if (($(tr).find("a.aslink").text()) == x || x == "all"){
            if (y) closebgctr($(tr));
            else openbgctr($(tr));
        }
    });
}

//Scroller to find bgc of interest
function findbgc(bgcId){
    openbgcrow(bgcId);
    $('html, body').animate(
    {
      scrollTop: $('#'+bgcId+'-svg').offset().top-120,
    },
    500,
    'linear');
}

aswindow = null;
function openaslink(x){
    if (aswindow) aswindow.close();
    aswindow = window.open(x.href, "bgcwindow");
}

function changegnselector(newvalue){
    $('#genetreeselector').val(newvalue);
    $('#genetreeselector').change();
}

$('#genetreeselector').on('change',function changegntree(){
    if (this.value.length) $('#gntreediv').html("<img src='trees/"+this.value+"' id='gntreeimg' width='800'>");
});

function copyToClipboard (text) {
    window.prompt ("Copy to clipboard: Ctrl+C, Enter", text);
}

bgctc=0
function toggleallbgcrow(e){
    e.preventDefault();
    openbgcrow("all",bgctc%2);
    bgctc+=1
}

function togglelog(){
    logwell = $("#logwell");
    if (logwell.css("display").toString()=="none"){
        logwell.css("display","inherit");
    }
    else{
        logwell.css("display","none");
    }
}

//Toggle column
$('a.toggle-vis').on( 'click', function (e) {
        e.preventDefault();
        // Get the column API object
        var column = summTable.column( $(this).attr('data-column') );
        // Toggle the visibility
        column.visible( ! column.visible() );
        $(".headertitles").tooltip();
    } );

//Toggle column kr table
$('a.toggle-vis2').on( 'click', function (e) {
        e.preventDefault();
        // Get the column API object
        var column = krTable.column( $(this).attr('data-column') );
        // Toggle the visibility
        column.visible( ! column.visible() );
        $(".headertitles").tooltip();
    } );

//Detect url params to go to bgc
var bgcIdparam = getUrlParameter("bgcid");
var tries = 0
var gotobgctimer;
function gotobgc(){
    tries+=1;
    console.log("trying to openbgc..."+tries);
    if (bgcIdparam && $('#'+bgcIdparam+'-link').length){
        findbgc(bgcIdparam);
        clearInterval(gotobgctimer);
    }
    else if (tries > 10){
        clearInterval(gotobgctimer);
    }
}
$(window).on("load",function(){
    if (bgcIdparam) gotobgctimer = setInterval(gotobgc,500);
});