export function size(ar){
  // get the dimensions of an array
  // source: https://stackoverflow.com/questions/49616639/how-can-i-export-all-functions-from-a-file-in-js
  var row_count = ar.length;
  var row_sizes = []
  for(var i=0;i<row_count;i++){
      row_sizes.push(ar[i].length)
  }
  return [row_count, Math.min.apply(null, row_sizes)]
}