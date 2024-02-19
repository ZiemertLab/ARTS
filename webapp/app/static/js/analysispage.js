// generates hashes used for linking multiselects and adding/removing from these
function hashCode(input) {
	var hash = 0;
	var j;
	if (input.length == 0) return hash;
	for (j = 0; j < input.length; j++) {
		c = input.charCodeAt(j);
		hash = ((hash<<5)-hash)+c;
		hash = hash & hash; // Convert to 32bit integer
	}
	return hash;
}

// check if number of seqs is over maximum
//Impose genome limits
function maxSeqs() {
if ($("#upfiles > option").length >= 200) {
            $("#seqbtn").addClass("disabled");
            $("#ncbibtn").addClass("disabled");
            $("#seqwarn").removeClass("hidden");
            return false;
        } else {
            $("#seqbtn").removeClass("disabled");
            $("#ncbibtn").removeClass("disabled");
            $("#seqwarn").addClass("hidden");
            return true;
        }
}

// prevent submitting with empty list
// check if any orgs in list
function validateForm() {
if ($('.selectablein').has('option').length>0) {
return "validated";
} else {
return "notvalidated";}
}
// disables submit button if list is empty
function submitToggle() {
    if ($('#upfiles').has('option').length>0) {
        $('#submitjob').removeClass('disabled');
    } else {
        $('#submitjob').addClass('disabled');
    }
}

function clearProgress() {
if ($("#uploadprog").attr("value") == $("#uploadprog").attr("max")) {
    $("#uploadprog").attr("value",0);
}
}

function uploadSuccess(data,textStatus,xhr) { //communication success
var data = JSON.parse(data);
console.log(data);
$("#seqbtn").removeClass("disabled");
$("#ncbibtn").removeClass("disabled");
$("#ncbiload").addClass("hidden");
$('progress').removeClass("hidden");
if (data["filename"] == false) {
    $("#errorwarning").removeClass("hidden");
    } else {
       if ($("#filesrc").val() == 'ncbi') {
        $("#uploadprog").attr("value",100);
        }
    for (var fileNumber = 0; fileNumber < data["filename"].length; fileNumber++) {
    $("#upfiles").prepend("<option value='"+data["filename"][fileNumber]+"' class='"+hashCode(data["filename"][fileNumber])+" picked'>"+data["name"][fileNumber]+"</option>");
    }
    maxSeqs();
    $("#uploadsuccess").removeClass("hidden");
    submitToggle();
    }
    clearProgress();
}
function uploadError(xhr,ajaxOptions,thrownError) { // communication error
console.log(thrownError);
$("#seqbtn").removeClass("disabled");
$("#ncbibtn").removeClass("disabled");
//unhide error message
$("#uploadwarning").removeClass("hidden");
$("#ncbiload").addClass("hidden");
clearProgress();
$("#uploadprog").attr("value",0);
}
function clearErrors(errorMessage) {
$(errorMessage).addClass("hidden");
}
function uploader(filesrc) {
if (maxSeqs()) {
$("#seqbtn").addClass("disabled");
$("#ncbibtn").addClass("disabled");
var uploadForm = {
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/upload',
        async: true,
        type: 'POST',
        //Form data
        data: new FormData($("#sequpload")[0]),
        cache: false,
        contentType: false,
        processData: false,
        success: uploadSuccess,
        error: uploadError
        /*xhr: function() {
            var myXhr = $.ajaxSettings.xhr();
            if (myXhr.upload) {
                // For handling the progress of the upload
                myXhr.upload.addEventListener('progress', function(e) {
                    if (e.lengthComputable) {
                        $('progress').attr({
                            value: e.loaded,
                            max: e.total,
                        });
                    }
                } , false);
            }
            return myXhr;
        }*/
};
//console.log(uploadForm.data);
if (filesrc == 'ncbi') {
    $('#ncbiload').removeClass("hidden");
    $('progress').addClass("hidden");
    $('progress').attr({
    value: 50,
    max: 100
    });
} else {
uploadForm.xhr = function() {
            var myXhr = $.ajaxSettings.xhr();
            if (myXhr.upload) {
                // For handling the progress of the upload
                myXhr.upload.addEventListener('progress', function(e) {
                    if (e.lengthComputable) {
                        $('progress').attr({
                            value: e.loaded,
                            max: e.total,
                        });
                    }
                } , false);
            }
            return myXhr;
        }
}
$.ajax(uploadForm);
}
}

function validateAcc() {
var userAcc = $("#ncbiacc1").val();
var allowedChars = /[A-Z0-9_]+$/i;
var wgs = /[A-Z]{4}0+(\.\d){0,}$/i;
var assembly = /^GC._/i;

console.log(wgs.test(userAcc));
console.log(assembly.test(userAcc));

if (userAcc && (userAcc.search(" ") == -1) && allowedChars.test(userAcc) && !(wgs.test(userAcc)) && !(assembly.test(userAcc))) { //check if userAcc is alphanumeric
    return true;
    } else {
    return false;
    }
}

function checkDuplicateAcc() {
var userAcc = $('#ncbiacc1').val().trim();
var duplicateFound;
$('#upfiles>option').each(function() {
    var plainName = $(this).text().split(".");
    var fileName = plainName[0];
    if (userAcc && userAcc.toUpperCase() == fileName.toUpperCase()) {
        duplicateFound=true;
    }
});
if (duplicateFound == true) {
    return false;
} else {
    return true;
}
}

function checkDuplicateFile() {
var listLoadedFiles = [];
$('#upfiles>option').each(function() {
    listLoadedFiles.push($(this).text());
})
var userFile = document.getElementById("seqfile1").files;
for (var counter = 0; counter < userFile.length; counter ++) {
    if ($.inArray(userFile[counter].name,listLoadedFiles) != -1) {
    return false;
    }
}
return true;
}

function uploadSequence(filesrc) {
clearErrors('#uploadwarning');
clearErrors('#errorwarning');
clearErrors('#uploadsuccess');
clearErrors('#accwarning');
clearErrors('#duplicatewarning');
clearErrors('.selectablewarn');
clearErrors('.flashes');
$("#filesrc").val(filesrc);
//checkDuplicateFile();
if (filesrc == 'ncbi' && !(validateAcc())) {
    $("#accwarning").removeClass("hidden");
} else if ((filesrc == 'ncbi' && !(checkDuplicateAcc())) || (filesrc == 'seqfile' && !(checkDuplicateFile()))) {
    $('#duplicatewarning').removeClass("hidden");
}
else {
uploader(filesrc);
}
}

function genusSuccess(data,textStatus,xhr) {
var counter;
for (counter=0;counter<data.length;counter++) {
    var genusInfo = data[counter];
    $('#genusselect').append("<option id='"+genusInfo.genusid+"' value='"+genusInfo.genusid+"'>"+genusInfo.genusname+"</option>");
}
}

function genusError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

function refGenus() {
$.ajax({
    contentType: 'application/json',
    url: '/results2/refgenus',
    async: true,
    cache: false,
    contentType: false,
    processData: false,
    success: genusSuccess,
    error: genusError

});
}

function removeAllSeqs() {
    var selectedValues3=$("#upfiles").val();
    var l;
    for (l = 0; l < selectedValues3.length;l++) {
    removeFromList('#upfiles', selectedValues3[l]);
    }
    maxSeqs();
    submitToggle();
}
function selectToggle() {
if (($('input[name="workflow"]:checked', '#sequpload').val()) == "2") {
    $('#genusform').removeClass("hidden");
    } else {
    $('#genusform').addClass("hidden");
    }
}

// check if accession number has been entered before upload button is available - currently deactivated

/*$('#ncbiacc1').on('change', function(){
    if ($('#ncbiacc1').val().length) {
        $('#ncbibtn').removeClass('disabled');
    } else {
        $('#ncbibtn').addClass('disabled');
    }
});*/

// check if file has been entered before upload button is available - currently deactivated due to compatibility issues with Safari

/*$('#seqfile1').on('change', function() {
    if ($('#seqfile1').val().length) {
        $('#seqbtn').removeClass('disabled');
    } else {
        $('#seqbtn').addClass('disabled');
    }
});*/

//Impose genome limits
function limitseqs(objid){
    var allowed = 200 - $("#upfiles > option").length;
    var allfiles = $("#"+objid)[0].files;
    var totalsize = 0;
    for (var i = 0; i < allfiles.length; ++i) {
        totalsize += allfiles.item(i).size;
        }
    //Check if seqs are less than total allowed and if request is under 100MB
    //if ($("#"+objid)[0].files.length > allowed || totalsize > 104857600){
    if ($("#"+objid)[0].files.length > allowed || totalsize > 1048576000000){
        $("#"+objid).val("");
        $("#uploadwarn").removeClass("hidden")
    }
    else{
        $("#uploadwarn").addClass("hidden")
    }
}

function toggleuploadbtn(btnid,x) {
    if (x.toString().length >= 1) {
        $("#"+btnid).prop("disabled",false);
    }
    else {
        $("#"+btnid).prop("disabled",true);
    }
}

$(document).ready(function(){
    $('[data-toggle="tooltip"]').tooltip();
});
$('[data-toggle="tooltip"]').on('click',function(){
    $(this).blur();
})
// prevent submitting form on hitting enter
$(document).on("keydown","input:not(button)", function(e) {
    if (e.keyCode == 13) {
        e.preventDefault();
    }
});

function validateForm(){
    return true;
}

// for submitting forms where multiselect is used as list of selected entries - all entries (under the cutoff) in the multiselect need to be selected first
function selectAndSend(formid) {
$(".picked").each(function() { this.selected=true; });
$(formid).submit();
}