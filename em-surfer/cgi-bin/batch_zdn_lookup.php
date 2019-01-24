<?php
include ('../header.php');
?>

<br>
<script type='text/javascript' src='../jsFunctions.js'> </script>
<script type='text/javascript' src='../ajaxFunctions-jsmol.js'> </script>
<script type='text/javascript' src='../jquery.all.js'> </script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.6.0/Chart.min.js"> </script>
<script type="text/javascript" src="../jsmol/JSmol.min.js"> </script>
<link rel="stylesheet" type="text/css" href="../jquery.autocomplete.css" />
<link rel="stylesheet" type="text/css" href="default.css" />
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">

<script>
var countDownDate = new Date().getTime() + 601000;
var x = setInterval(function() {
    var now = new Date().getTime();
    var distance = countDownDate - now;
    var days = Math.floor(distance / (1000 * 60 * 60 * 24));
    var hours = Math.floor((distance % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60));
    var minutes = Math.floor((distance % (1000 * 60 * 60)) / (1000 * 60));
    var seconds = Math.floor((distance % (1000 * 60)) / 1000);
    
    if (seconds >= 10) document.getElementById("timer").innerHTML = "The zip file is avaliable for: " + minutes + ": " + seconds;
    else document.getElementById("timer").innerHTML = "The zip file is avaliable for: " + minutes + ": 0" + seconds;
    
    if (distance < 0) {
        clearInterval(x);
        document.getElementById("download_button").innerHTML = "";
        document.getElementById("timer").innerHTML = "The zip file is expired, please perform the search again. ";
    }
}, 1000);
</script>

<?php

function generateRandomString($length) {
    $characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
    $charactersLength = strlen($characters);
    $randomString = '';
    for ($i = 0; $i < $length; $i++) {
        $randomString .= $characters[rand(0, $charactersLength - 1)];
    }
    return $randomString;
}

$query_pdb_id = $_POST['listarea'];
$query_pdb_file = $_FILES['listfile']['tmp_name'];
$volume_filter = $_POST['volumeFilter'];
$num = $_POST['num'];
$representation = $_REQUEST['rep'];

$inv_dir = '../db/2_contour/inv';
if ($representation === "recommend") {
    $inv_dir = '../db/recommend_contour/inv';
}
if ($representation === "recommendonethird") {
    $inv_dir = '../db/1_contour/inv';
}
if ($representation === "recommendtwothird") {
    $inv_dir = '../db/2_contour/inv';
}
if ($representation === "recommendonethirdtwothird") {
    $inv_dir = '../db/1_2_contour/inv';
}

$zernike_length_end_pos;
$zernike_descriptors;
$zernike_count;
$current_inv_directory;

$query_pdb_file = $_FILES['listfile']['tmp_name'];
if ($query_pdb_file) {
    $unique_filename = substr($query_pdb_file, 5);
    #move
    $target = "upload/" . $unique_filename;
    move_uploaded_file($query_pdb_file, $target);

    # read the content of the uploaded file
    $fh = fopen($target, 'r');
    $thefiledata = fread($fh, filesize($target));
    fclose($fh);

    $query = $thefiledata;
} else {
    $query = $query_pdb_id;
}

$current_inv_directory = $inv_dir;
$query_file = "$current_inv_directory/$query.inv";
//echo $query_file;

$pattern = "/:|\s+|\n/";
$query = preg_replace($pattern, ",", $query);

$query_list = explode(",", $query); 
$query_item_existance = array(); 
$zdn_list = array(); 
$structure_list = array(); 
$output_file_list = array(); 
$output_file_string = ""; 

for ($i = 0; $i < count($query_list); $i++) {
    $query_file = "$current_inv_directory/$query_list[$i].inv"; 
    file_exists($query_file) ? array_push($query_item_existance, true) : array_push($query_item_existance, false); 
}

for ($i = 0; $i < count($query_list); $i++) {

    array_push($structure_list, "$query_list[$i]"); 
    if (!$query_item_existance[$i]) {
        array_push($zdn_list, ""); 
        continue; 
    } 

    if (in_array($query_list[$i] . '.txt', $output_file_list)) continue; 

    $zernike_count = 0;
    $descriptor;
    $textarea_formatted_descriptors = "";
    $textarea_formatted_output = ""; 

    $query_file = "$current_inv_directory/$query_list[$i].inv"; 
    $cat_inv = file_get_contents($query_file);
    $zernike_length_end_pos = strpos($cat_inv, "\n");
    $zernike_descriptors = substr($cat_inv, $zernike_length_end_pos);

    $tmparray = preg_split("/\n/", $zernike_descriptors);
    for ($j = 0; $j < count($tmparray) - 1; ++ $j) {
        $descriptor = $tmparray[$j];
        if ($zernike_count > 0) {
            $textarea_formatted_descriptors = "$textarea_formatted_descriptors($zernike_count) $descriptor" . " <br />";
            $textarea_formatted_output = "$textarea_formatted_output($zernike_count) $descriptor" . " \r\n";
        }
        $zernike_count++;
    }

    $file_name = $query_list[$i] . '.txt'; 
    array_push($output_file_list, $file_name); 
    $output_file_string = $output_file_string . $file_name . " "; 

    file_put_contents('queryresult/' . $file_name, $textarea_formatted_output);  
    array_push($zdn_list, $textarea_formatted_descriptors); 
}

$zip_file_name = 'result_' . generateRandomString(10) . '.zip'; 
shell_exec('cd queryresult && zip '. $zip_file_name . ' ' . $output_file_string); 

?>

<div class="jumbotron" id="zdn_batch_download">
  <h2>Your 3DZD result is ready!</h2>
  <p class="lead">Click the button below to download the zip file which contains all 3DZD info below.</p>
  <hr style="color: black;">
  <p id="timer">The zip file is avaliable for: 10:00</p>
  <div id="download_button">
    <a href="queryresult/<?php echo $zip_file_name ?>" id="zdn_batch_download_button" class="btn btn-dark btn-md" href="#" role="button">Download The 3DZD Compressed Zip Package</a>
  </div>
</div>

<div class="strike">
    <span id="strike_text">INDIVIDUAL RESULT BELOW</span>
</div>

<table class="table table-bordered table-sm" width="760px" align="center" border="1" id="zdn_batch_result">
    <?php 
        foreach($zdn_list as $index => $zdn_descriptor):
        if ($query_item_existance[$index]) {
    ?>
    <tr>
        <td id="zdn_batch_exist">
            <?php echo "$structure_list[$index]"; ?>
        </td>
        <td>
            3D Zernlike Descriptors
        </td>
    </tr>
    <tr>
        <td colspan="6">
            <div id="zdn_batch_content">
                <?php echo "$zdn_descriptor"; ?>
            </div>
        </td>
    </tr>
    <?php } else { ?>
    <tr>
        <td id="zdn_batch_not_exist">
            <?php echo "$structure_list[$index]"; ?>
        </td>
        <td>
            This structure does not exist 
        </td>
    </tr>
<?php } endforeach; ?>
</table>

<?php 
    include ('../tailer.php');  
    shell_exec('./discard_zdn_batch_zip.sh ' . $zip_file_name . ' >/dev/null 2>/dev/null & '); 
?>
