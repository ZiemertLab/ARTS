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
var summGenomes = 0;
var krTable = 0;
var dupTable = 0;
var bgcTable = 0;
var corecircle = 0;
var summchart = 0;
var statustimer = 0;
var genetreetimer = 0;
var finished = 0;
var all_regions = 0;
var filetojobid = {};

var sortgeneselector = function() {
    var optarr = $("#genetreeselector option");
    optarr.sort(function(a,b) {
        if (a.text.toLowerCase() > b.text.toLowerCase()) return 1;
        if (a.text.toLowerCase() < b.text.toLowerCase()) return -1;
        return 0; })
    $("#genetreeselector").empty().append(optarr);
}

function updatestatus(data){
    if (data.state){
        $('#statetxt').html(data.state);
    }
    pwidth = "1%"
    if (data.state.toLowerCase() == "running bgc networking"){
            pwidth = "80%"
    }
    if (data.state.toLowerCase() == "done"){
            $('#jobprogbar').addClass("progress-bar-success");
            $('#jobprogbar').removeClass("progress-bar-striped");
            $('#jobprogbar').removeClass("active");
            pwidth = "100%";
            //Update bigscape section
            $('#bigscapelink').html('<a class="btn btn-primary" style="width:400px" href="bigscape/" target="_blank">Open bigscape overview</a><br>');
            clearInterval(statustimer);
        }
    $('#jobprogbar').width(pwidth);

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
}
var summGenomes;
//DATATABLES
$(document).ready(function(){
    summGenomes = $('#genomeTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "genometable"
            },
//        dom: "t<'row'<'col-sm-12'f>>" +
//            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
//        pagingType: "simple",
        paging: false,
        searching: false,
        lengthChange: false,
        info: false,
        columns: [
            {"data":0},
            {"data":1, "visible":false,"defaultContent":""},
            {"data":2},
            {"data":3},
            {"data":4},
            {"data":5},
            {"data":6},
            {"data":7},
            {"data":8},
            {"data":9},
            {"data":10}
        ],
        createdRow: function ( row, data, index ) {
            $('td',row).eq(0).html("<a href='/results/"+data[1]+"' target='newartswindow'>"+data[0]+"</a>");
            filetojobid[data[0]]=data[1];
        },
        autoWidth: false
    });
    krTable = $('#krTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "multikrtab"
            },
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        columns: [
//            {
//                "className":      'details-control',
//                "orderable":      false,
//                "data":           null,
//                "defaultContent": '',
//                "width":"25px"
//            },
            {"data":0},
            {"data":1},
            {"data":2}
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
            "url": "multiduptab"
            },
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        columns: [
//            {
//                "className":      'details-control',
//                "orderable":      false,
//                "data":           null,
//                "defaultContent": '',
//                "width":"25px"
//            },
            {"data":0},
            {"data":1},
            {"data":2}
        ],
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
        autoWidth: false
    });
    reloadbgctimer = 0;
    bgcTable = $('#bgcTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "multibgctab"
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
            {"data":0, "width":"25px"},
            {"data":1, "width":"25px"},
            {"data":2, "render": function ( data, type, row ){
                    newdata = [];
                    for (i=0; i<data.length; i++){
                        fname = row[6][i].split("/").pop();
                        bgc = row[11][i].replace("cluster-","r").replace("_","c");
                        if(filetojobid.hasOwnProperty(fname))
                            newdata.push("<a href='/results/"+filetojobid[fname]+"/report?bgcid="+ bgc +"' target='newartswindow'>"+data[i]+"</a>");
                        else{
                            newdata.push(data[i]);
                            clearTimeout(reloadbgctimer);
                            reloadbgctimer = setTimeout(function(){ bgcTable.ajax.reload(); console.log("Reloaded BGC table");},2000);
                        }
                    }
                    return newdata.join("<br>");
                }},
            {"data":3, "render": function ( data, type, row ){ return data.join("<br>")}},
            {"data":4, "visible":false},
            {"data":7, "width":"25px"},
            {"data":8, "width":"25px"},
            {"data":9, "visible":false}
        ],
        order: [[2,'desc'],[6,'desc']],
        createdRow: function ( row, data, index ) {
            if (data[7]>0){
                $('td',row).eq(5).addClass("summtab_yes");
            }
            if (data[8]>0){
                $('td',row).eq(6).addClass("summtab_yes");
            }
        },
        aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
        autoWidth: false
    });
    summTable = $('#summaryTable').DataTable({
        ajax: {
            "contentType": "application/json",
            "url": "multicoretab"
            },
        order: [5,'desc'],
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
            {"data":5},
            {"data":6},
            {"data":7, visible:false},
            {"data":8, visible:false}
        ],
        dom: "t<'row'<'col-sm-12'f>>" +
            "<'row'<'col-sm-4'l><'col-sm-4'i><'col-sm-4'p>>B",
        pagingType: "simple",
        createdRow: function ( row, data, index ) {
            if (data[3] != "0.0"){
                $('td',row).eq(4).addClass("summtab_yes");
            }
            if (data[4] != "0.0"){
                $('td',row).eq(5).addClass("summtab_yes");
            }
            if (data[5] != "0.0"){
                $('td',row).eq(6).addClass("summtab_yes");
            }
            if (data[6] != "0.0"){
                $('td',row).eq(7).addClass("summtab_yes");
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
    $.getJSON('multistatus',updatestatus);
    clearInterval(statustimer);
    statustimer = setInterval(function() {
        $.getJSON('multistatus',updatestatus)
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


/* Formatting function for row details */
function formatbgcrow ( d ) {
    childstring = '<div class="bgctabChild"><h5>BGCs</h5><table cellpadding="5" class="table table-bordered table-striped">';
    childstring += "<tr><th>Cluster</th><th>Source</th><th>Core Hits</th><th>Res Hits</th></tr>";

    clusters = d[2];
    clustIDS = d[11];
    sources = d[6];
    corehits = d[9];
    reshits = d[10];
//    clusters = JSON.parse(d[2].replace(/u'/g,'"').replace(/'/g,'"'));
//    clustIDS = JSON.parse(d[11].replace(/u'/g,'"').replace(/'/g,'"'));
//    sources = JSON.parse(d[6].replace(/u'/g,'"').replace(/'/g,'"'));
//    corehits = JSON.parse(d[9].replace(/u'/g,'"').replace(/'/g,'"'));
//    reshits = JSON.parse(d[10].replace(/u'/g,'"').replace(/'/g,'"'));

    hitstoclusts = {}

    for (i=0; i < clusters.length; i++){
        fname = sources[i].split('/').pop();
        sourcestring = fname;
        if (fname in filetojobid){
            sourcestring = "<a href='/results/"+filetojobid[fname]+"' target='newartswindow'>"+fname+"</a>";
        }
        childstring += "<tr><td>"+"<a href='/results/"+filetojobid[fname]+"/report?bgcid="+clustIDS[i].replace(/cluster-/g,"r").replace(/_/g,"c")+"' target='newartswindow'>"+clusters[i]+"</a></td>";
        childstring += "<td>"+sourcestring+"</td><td>"+corehits[i].join("<br>")+"</td><td>"+reshits[i].join("<br>")+"</td></tr>";
        hits = corehits[i].concat(reshits[i]);
        for (j=0; j < hits.length; j++){
            if (!(hits[j] in hitstoclusts)){
                hitstoclusts[hits[j]]=[]
            }
            hitstoclusts[hits[j]].push("<a href='/results/"+filetojobid[fname]+"/report?bgcid="+clustIDS[i].replace(/cluster-/g,"r").replace(/_/g,"c")+"' target='newartswindow'>"+clusters[i]+"</a>");
        }
    }

    childstring += '</table><br><h5>Hits in family</h5><table cellpadding="5" class="table table-bordered table-striped"><tr><th>Model Hit</th><th>Clusters</th><th>Family presence</th></tr>';

    for (var key in hitstoclusts) {
        if (hitstoclusts.hasOwnProperty(key)) {
            childstring += '<tr><td>'+key+'</td><td>'+hitstoclusts[key].join("<br>")+'</td><td>'+hitstoclusts[key].length+" / "+clusters.length+'</td></tr>'
        }
    }

    childstring += "</table></div>"

    return childstring;
}

var sourcetobgc = null;

function formatchildrow ( d ) {
    //Create lookup table if not present
    if (sourcetobgc == null){
        bgcjson = bgcTable.ajax.json();
        sourcetobgc = {};
        for (i=0; i<bgcjson["data"].length; i++){
            b = bgcjson["data"][i];

            clusters = b[2];
            clustIDS = b[11];
            gcfs = b[6];
            corehits = b[9];
            for (j=0; j<gcfs.length; j++){
                gcf = gcfs[j].split('/').pop();
                if (sourcetobgc.hasOwnProperty(gcf) == false){
                    sourcetobgc[gcf] = {};
                }
                //index by core hits
                for (c=0; c<corehits[j].length; c++){
                    core = corehits[j][c];
                    if (sourcetobgc[gcf].hasOwnProperty(core) == false){
                        sourcetobgc[gcf][core] = [];
                    }
                    sourcetobgc[gcf][core].push("<a href='/results/"+filetojobid[gcf]+"/report?bgcid="+clustIDS[j].replace(/cluster-/g,"r").replace(/_/g,"c")+"' target='newartswindow'>"+clusters[j]+"</a> (GCF: "+b[0]+")");
                }
            }
        }
    }
    childstring = '<table cellpadding="5" class="table summtabChild">';
    childstring += '<tr><th>Source</th><th>BGC(s)</th><th>Seq(s) (start,end,strand)</th><th>Dup</th><th>BGC</th><th>HGT</th><th>Res</th></tr>';
    for (i = 0; i < d[8].length; i++){
        sourcestring = d[8][i];
        //Find corresponding BGCs
        bgcs = [""];
        if (sourcetobgc.hasOwnProperty(sourcestring) && sourcetobgc[sourcestring].hasOwnProperty(d[0])) bgcs = sourcetobgc[sourcestring][d[0]];

        // Find source in other hits lists
        isdup = "";
        if(d[9].join(",").search(sourcestring)!= -1) isdup = "<span class='glyphicon glyphicon-ok'></span>";
        isbgc = "";
        if(d[10].join(",").search(sourcestring)!= -1) isbgc = "<span class='glyphicon glyphicon-ok'></span>";
        ishgt = "";
        if(d[11].join(",").search(sourcestring)!= -1) ishgt = "<span class='glyphicon glyphicon-ok'></span>";
        isres = "";
        if(d[12].join(",").search(sourcestring)!= -1) isres = "<span class='glyphicon glyphicon-ok'></span>";

        if (d[8][i] in filetojobid){
            sourcestring = "<a href='/results/"+filetojobid[d[8][i]]+"' target='newartswindow'>"+d[8][i]+"</a>"
        }
        childstring += "<tr><td>"+sourcestring+"</td><td>"+bgcs.join("<br>")+"</td><td>"+d[13][i].join("<br>")+"</td><td>"+isdup+"</td><td>"+isbgc+"</td><td>"+ishgt+"</td><td>"+isres+"</td></tr>";
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
    rslt = formatbgcrow(d);
    row.child( rslt ).show();
    tr.addClass('shown');
//    if (clust) svgene.drawClusters(d[0].replace("cluster-","r").replace("_","c")+"-svg", [rslt[1]], 20, 800);
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