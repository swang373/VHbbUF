<html>
<body>

<form action="b.php" method="post">
<input type="submit" name="submit" value="REFRESH"/>
</form>

<?php
$var = $_POST['var'];
$xtitle   = $_POST['xtitle'];
$cuts   = $_POST['cuts'];
$xlow   = $_POST['xlow'];
$xhigh   = $_POST['xhigh'];
$binsize   = $_POST['binsize'];
$log = 1;

if(isset($_POST['formWheelchair']) && $_POST['formWheelchair'] == 'Yes')
{
 $log = 1;
}
else 
{
 $log = 0;
}

$formatted = sprintf("echo %s %s >> PISS", $fname, $age);
$proc=proc_open("tail PISS",
  array(
    array("pipe","r"),
    array("pipe","w"),
    array("pipe","w")
  ),
  $pipes);
print stream_get_contents($pipes[1]);
?>
</body>
</html>


