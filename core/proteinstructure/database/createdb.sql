SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL';

DROP SCHEMA IF EXISTS `KiharaDB` ;
CREATE SCHEMA IF NOT EXISTS `KiharaDB` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci ;
USE `KiharaDB` ;

-- -----------------------------------------------------
-- Table `KiharaDB`.`StructureType`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `KiharaDB`.`StructureType` ;

CREATE  TABLE IF NOT EXISTS `KiharaDB`.`StructureType` (
  `StructureTypeID` TINYINT NOT NULL ,
  `Description` VARCHAR(20) NOT NULL ,
  PRIMARY KEY (`StructureTypeID`) ,
  UNIQUE INDEX `UNIQUE_Description` (`Description` ASC) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `KiharaDB`.`Structure`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `KiharaDB`.`Structure` ;

CREATE  TABLE IF NOT EXISTS `KiharaDB`.`Structure` (
  `StructureID` BIGINT NOT NULL ,
  `Description` VARCHAR(20) NULL ,
  `StructureFilePath` VARCHAR(1000) NULL ,
  `ImageFilePath` VARCHAR(1000) NULL ,
  `AnimationFilePath` VARCHAR(1000) NULL ,
  `ZernikeDescriptors` VARCHAR(1000) NULL ,
  `ParentStructureID` BIGINT NULL ,
  `StructureTypeID` TINYINT NOT NULL ,
  PRIMARY KEY (`StructureID`) ,
  INDEX `INDEX_Structure_Structure_ParentStructureID_StructureID` (`ParentStructureID` ASC) ,
  INDEX `INDEX_Structure_StructureType_StructureTypeID_StructureTypeID` (`StructureTypeID` ASC) ,
  CONSTRAINT `FK_Structure_Structure_ParentStructureID_StructureID`
    FOREIGN KEY (`ParentStructureID` )
    REFERENCES `KiharaDB`.`Structure` (`StructureID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `FK_Structure_StructureType_StructureTypeID_StructureTypeID`
    FOREIGN KEY (`StructureTypeID` )
    REFERENCES `KiharaDB`.`StructureType` (`StructureTypeID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `KiharaDB`.`MultipleParentStructures`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `KiharaDB`.`MultipleParentStructures` ;

CREATE  TABLE IF NOT EXISTS `KiharaDB`.`MultipleParentStructures` (
  `StructureID` BIGINT NOT NULL ,
  `ParentStructureID` BIGINT NOT NULL ,
  PRIMARY KEY (`StructureID`, `ParentStructureID`) ,
  INDEX `INDEX_MultipleParentStructures_Structure_StructureID_StructureID` (`StructureID` ASC) ,
  INDEX `INDEX_MultipleParentStructures_Structure_ParentID_StructureID` (`ParentStructureID` ASC) ,
  CONSTRAINT `FK_MultipleParentStructures_Structure_StructureID_StructureID`
    FOREIGN KEY (`StructureID` )
    REFERENCES `KiharaDB`.`Structure` (`StructureID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `FK_MultipleParentStructures_Structure_ParentID_StructureID`
    FOREIGN KEY (`ParentStructureID` )
    REFERENCES `KiharaDB`.`Structure` (`StructureID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `KiharaDB`.`CATH`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `KiharaDB`.`CATH` ;

CREATE  TABLE IF NOT EXISTS `KiharaDB`.`CATH` (
  `CATHID` INT NOT NULL ,
  `Class` INT NOT NULL ,
  `Architecture` INT NOT NULL ,
  `Topology` INT NOT NULL ,
  `Homologous` INT NOT NULL ,
  PRIMARY KEY (`CATHID`) ,
  UNIQUE INDEX `UNIQUE_Class_Architecture_Topology_Homologous` (`Class` ASC, `Architecture` ASC, `Topology` ASC, `Homologous` ASC) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `KiharaDB`.`StructureCATHCodes`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `KiharaDB`.`StructureCATHCodes` ;

CREATE  TABLE IF NOT EXISTS `KiharaDB`.`StructureCATHCodes` (
  `StructureID` BIGINT NOT NULL ,
  `CATHID` INT NOT NULL ,
  PRIMARY KEY (`StructureID`, `CATHID`) ,
  INDEX `INDEX_StructureCATHCodes_Structure_StructureID_StructureID` (`StructureID` ASC) ,
  INDEX `INDEX_StructureCATHCodes_CATH_CATHID_CATHID` (`CATHID` ASC) ,
  CONSTRAINT `FK_StructureCATHCodes_Structure_StructureID_StructureID`
    FOREIGN KEY (`StructureID` )
    REFERENCES `KiharaDB`.`Structure` (`StructureID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `FK_StructureCATHCodes_CATH_CATHID_CATHID`
    FOREIGN KEY (`CATHID` )
    REFERENCES `KiharaDB`.`CATH` (`CATHID` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


;

CREATE USER `webuser` IDENTIFIED BY 'webuser';

grant SELECT on TABLE `KiharaDB`.`StructureCATHCodes` to webuser;
grant SELECT on TABLE `KiharaDB`.`CATH` to webuser;
grant SELECT on TABLE `KiharaDB`.`StructureType` to webuser;
grant SELECT on TABLE `KiharaDB`.`MultipleParentStructures` to webuser;
grant SELECT on TABLE `KiharaDB`.`Structure` to webuser;
;

CREATE USER `autoupdateuser` IDENTIFIED BY 'autoupdateuser';

grant INSERT on TABLE `KiharaDB`.`StructureCATHCodes` to autoupdateuser;
grant UPDATE on TABLE `KiharaDB`.`StructureCATHCodes` to autoupdateuser;
grant DELETE on TABLE `KiharaDB`.`StructureCATHCodes` to autoupdateuser;
grant UPDATE on TABLE `KiharaDB`.`CATH` to autoupdateuser;
grant INSERT on TABLE `KiharaDB`.`CATH` to autoupdateuser;
grant DELETE on TABLE `KiharaDB`.`CATH` to autoupdateuser;
grant UPDATE on TABLE `KiharaDB`.`StructureType` to autoupdateuser;
grant INSERT on TABLE `KiharaDB`.`StructureType` to autoupdateuser;
grant DELETE on TABLE `KiharaDB`.`StructureType` to autoupdateuser;
grant DELETE on TABLE `KiharaDB`.`MultipleParentStructures` to autoupdateuser;
grant INSERT on TABLE `KiharaDB`.`MultipleParentStructures` to autoupdateuser;
grant UPDATE on TABLE `KiharaDB`.`MultipleParentStructures` to autoupdateuser;
grant DELETE on TABLE `KiharaDB`.`Structure` to autoupdateuser;
grant INSERT on TABLE `KiharaDB`.`Structure` to autoupdateuser;
grant UPDATE on TABLE `KiharaDB`.`Structure` to autoupdateuser;
grant SELECT on TABLE `KiharaDB`.`StructureCATHCodes` to autoupdateuser;
grant SELECT on TABLE `KiharaDB`.`CATH` to autoupdateuser;
grant SELECT on TABLE `KiharaDB`.`StructureType` to autoupdateuser;
grant SELECT on TABLE `KiharaDB`.`MultipleParentStructures` to autoupdateuser;
grant SELECT on TABLE `KiharaDB`.`Structure` to autoupdateuser;

SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
