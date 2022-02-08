var data = source.data;
var filetext = 'contig\tcluster\n';
for (var i = 0; i < data.contig.length; i++) {
    var currRow = [data.contig[i].toString(),
                   data.cluster[i].toString(),
                   ];

    var joined = currRow.join('\t');
    filetext = filetext.concat(joined + '\n');
}

var filename = 'contigs.list';
var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
} else {
    var link = document.createElement("a");
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob);
    link.download = filename
    link.target = "_blank";
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'))
}
