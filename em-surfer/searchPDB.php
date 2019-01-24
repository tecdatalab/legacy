<?php

if (isset($_GET['q']) && $_GET['q'] != '') {
	//Add slashes to any quotes to avoid SQL problems.
	$search = addslashes($_GET['q']);
	$result=`ls db/2_contour/inv`;
	$inv_id_array = preg_split("/\n/", $result);
	for ($i = 0; $i < count($inv_id_array) - 2; ++ $i) {
		$inv_id_array[$i] = substr($inv_id_array[$i], 0, strrpos($inv_id_array[$i] ,".inv"));
	}
	$k=0;
	for ($j = 0; $j < count($inv_id_array) - 2; ++ $j) {
		$position=strrpos($inv_id_array[$j] ,$search);
		if ($position=== false) {
		} else {
			if($position == 0) {
    				echo $inv_id_array[$j]. "\n";
    				++$k;
    				if($k > 10){
    					break;
    				}
    			}
		}
	}
}
?>
